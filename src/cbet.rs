use crate::beam::*;
use crate::mesh::*;
use crate::consts;
use std::thread;

/// Do the CBET calculation! Populates the i_b fields of each crossing
pub fn cbet(mesh: &Mesh, beams: &mut [Beam], nthreads: usize) {
    consts::init_consts(); // inits CS const
    let mut currmax = consts::MAX_INCR;
    create_raystore(beams, mesh.nx, mesh.nz);

    for i in 1..=500 {
        get_cbet_gain(mesh, beams, nthreads);
        let updateconv = update_intensities(beams, 0.0, currmax);
        if updateconv <= consts::CONVERGE {
            break;
        }

        let currmaxa = consts::MAX_INCR*f64::powi(consts::CBETCONVERGENCE, i);
        let currmaxb = consts::CBETCONVERGENCE*updateconv;
        currmax = f64::min(currmaxa, currmaxb);
    }
    println!("Finished cbet, abt to post");
    post(mesh, beams);
}

/// Updates the i_b values of each crossing, one last time! Yipee!
fn post(mesh: &Mesh, beams: &mut [Beam]) {
    let w0 = 2.0*std::f64::consts::PI*consts::C_SPEED/consts::LAMBDA;
    let norm_factor_const = f64::sqrt(8.0*std::f64::consts::PI/consts::C_SPEED) * consts::ESTAT / (consts::ME_G*consts::C_SPEED*w0) * f64::sqrt(1e14*1e7);

    beams.iter().for_each(|beam| {
        beam.rays.iter().for_each(|ray| {
            ray.crossings.iter().enumerate().for_each(|(crossingnum, crossing)| {
                let mut crossing = crossing.lock().unwrap();
                //let crossing_next = ray.crossings[usize::min(crossingnum+1, ray.crossings.len()-1)].lock().unwrap();
                let next_area_ratio = if crossingnum+1 == ray.crossings.len() {
                    crossing.area_ratio
                } else {
                    let crossing_next = ray.crossings[crossingnum+1].lock().unwrap();
                    crossing_next.area_ratio
                };
                // some of this code looks copied from the cbet math code (TODO)
                let ix = crossing.boxesx;
                let iz = crossing.boxesz;

                let area_avg = (crossing.area_ratio+next_area_ratio)/2.0;
                let ne_over_nc = mesh.get(ix, iz).eden/consts::NCRIT;
                let ne_over_nc_corrected = f64::min(ne_over_nc, 1.0);
                let ne_term = f64::sqrt(1.0-ne_over_nc_corrected);
                let epsilon_eff = ne_term*ne_term;
                let interaction_mult = 1.0/(area_avg*ne_term)*1.0/f64::sqrt(epsilon_eff);
                let norm_factor = norm_factor_const * f64::sqrt(interaction_mult) * epsilon_eff.powf(0.25);
                let intensity_new = f64::sqrt(crossing.i_b) * norm_factor;

                // now we have mutable borrow
                crossing.i_b = intensity_new;
            });
        });
    });
}

/// Updates the i_b values of each crossing. Returns the new value of convMax (which is
/// updateconv in cbet fn.). The variable names are different because that's how it is in the
/// c++ implementation.
fn update_intensities(beams: &mut [Beam], conv_max: f64, curr_max: f64) -> f64 {
    let mut curr_conv_max = conv_max;
    beams.iter().for_each(|beam| {
        beam.rays.iter().for_each(|ray| {
            let i0 = {
                let first_crossing = ray.crossings[0].lock().unwrap();
                first_crossing.i_b
            };
            let mut mult_acc = 1.0;
            ray.crossings.iter().for_each(|crossing| {
                let mut crossing = crossing.lock().unwrap();
                let (new_intensity, new_conv_max) = limit_energy(&crossing, mult_acc, i0, curr_max, curr_conv_max);
                curr_conv_max = new_conv_max;
                mult_acc *= crossing.w_mult;
                crossing.i_b = new_intensity;
            });
        });
    });
    curr_conv_max
}

/// limitEnergy fn. from cpp, returns new value of updateConv/convMax/maxChange, whichever you
/// choose to call it haha
///
/// Again, c++ variable names copied, some are different from what they're named by the
/// function that calls them.
fn limit_energy(crossing: &Crossing, multiplier_acc: f64, i0: f64, curr_max: f64, max_change: f64) -> (f64, f64) {
    // copied comments from c++ impl!
    let i_prev = crossing.i_b; // previous value in cell
    let mut i_curr = i0*multiplier_acc; // "current" value, unclamped update
    // the fractional change in energy from imposing the update as is
    let fractional_change = f64::abs(i_curr-i_prev)/i_prev;
    // update the convergence check variable
    let new_max_change = f64::max(fractional_change, max_change);
    // if the fractional change is too large, clamp the value
    if fractional_change > curr_max {
        let sign = if i_curr - i_prev > 0.0 { 1.0 } else { -1.0 };
        let correction = 1.0 + curr_max*sign;
        i_curr = i_prev*correction;
    }
    (i_curr, new_max_change)
}

/// Get CBET gain. Modifies the w_mult property of each crossing.
fn get_cbet_gain(mesh: &Mesh, beams: &mut [Beam], nthreads: usize) {
    beams.iter().enumerate().for_each(|(beamnum, beam)| {
        thread::scope(|s| {
            for chunk in beam.rays.chunks(beam.rays.len() / nthreads) {
                s.spawn(|| {
                    chunk.iter().for_each(|ray| {
                        ray.crossings.iter().for_each(|crossing| {
                            // create a clone of the crossing
                            let crossing_clone = {
                                let crossing_mguard = crossing.lock().unwrap();
                                crossing_mguard.clone()
                            };

                            let other_beam_cbet_incr = |other_beam: &Beam| {
                                // raycross/raycross_next are also clones
                                let (raycross, raycross_next) = {
                                    let (has_crossing, (other_rayind, raycross_ind)) =
                                        other_beam.raystore[crossing_clone.boxesx*mesh.nz+crossing_clone.boxesz];
                                    if !has_crossing {
                                        return 0.0;
                                    }
                                    let other_ray_crossings = &other_beam.rays[other_rayind].crossings;
                                    let raycross = other_ray_crossings[raycross_ind].lock().unwrap();
                                    let raycross_next = if raycross_ind+1 == other_ray_crossings.len() {
                                        raycross.clone()
                                    } else {
                                        let rcnext_mguard = other_ray_crossings[raycross_ind+1].lock().unwrap();
                                        rcnext_mguard.clone()
                                    };
                                    (raycross.clone(), raycross_next)
                                };

                                get_cbet_increment(mesh, &crossing_clone, &raycross, &raycross_next)
                            };
                            let cbet_sum = beams.iter().enumerate().map(|(otherbeamnum, other_beam)| {
                                if beamnum == otherbeamnum {
                                    return 0.0;
                                }
                                other_beam_cbet_incr(other_beam)
                            }).sum::<f64>();

                            let mut crossing = crossing.lock().unwrap();
                            crossing.w_mult = f64::exp(-1.0*cbet_sum);
                        });
                    });
                });
            }
        });
    });
}

/// This is where all of the math lives. Lines 162 through 209 of cbet.cpp are translated here.
/// Only thing that is changed is variable names are changed to snake case cuz I don't want
/// cargo to yell at me :(
fn get_cbet_increment(mesh: &Mesh, crossing: &Crossing, raycross: &Crossing, raycross_next: &Crossing) -> f64 {
    let mesh_pt = mesh.get(crossing.boxesx, crossing.boxesz);
    // I'm copying comments over from the c++ implementation. They're a little verbose
    // (even for me). Have fun reading!

    // INTERACTION MULTIPLIER: find ray_o's interaction multiplier
    // get the avereage area of the ray across the zone
    let area_avg = (raycross.area_ratio+raycross_next.area_ratio)/2.0;
    // NOTE: This neOverNc value can be taken straight from the grid
    let ne_over_nc = mesh_pt.eden/consts::NCRIT;
    // clamp the maximum neOverNc value to
    let ne_over_nc_corrected = f64::min(ne_over_nc, 1.0);
    // term used multiple times in calculations
    let ne_term = f64::sqrt(1.0-ne_over_nc_corrected);
    let epsilon_eff = ne_term*ne_term;
    let interaction_mult = 1.0/(area_avg*ne_term)*1.0/f64::sqrt(epsilon_eff);

    // eta Routine
    let kx_seed = crossing.dkx;
    let kz_seed = crossing.dkz;
    let kx_pump = raycross.dkx;
    let kz_pump = raycross.dkz;

    let machx = mesh_pt.machnum; // the mach number of the plasma velocity
    let machz = 0.0; // TODO: Implement multidimensional plasma velocity

    // assuming uniform wavelength/frequency
    let omega1 = consts::OMEGA;
    let omega2 = consts::OMEGA;
    // find ion acoustic wave vector, difference of ray trajectories scaled to
    // omega/(sqrt(1-neOverNc)*c)

    let iaw_vector = [
        (omega1*kx_seed - omega2*kx_pump)*f64::sqrt(1.0-ne_over_nc)/consts::C_SPEED,
        (omega1*kz_seed - omega2*kz_pump)*f64::sqrt(1.0-ne_over_nc)/consts::C_SPEED,
    ];
    let k_iaw = f64::sqrt(iaw_vector[0]*iaw_vector[0] + iaw_vector[1]*iaw_vector[1]); // magnitude of iaw vector
    let eta_numerator;
    let eta_denominator;
    // WILLOW SAYS: have to use CS which is a mutable static (global var) because it is
    // initialized at runtime (uses non-constant sqrt fn)
    unsafe {
        eta_numerator = omega1-omega2 - (iaw_vector[0]*machx + iaw_vector[1]*machz)*consts::CS;
        eta_denominator = k_iaw*consts::CS;
    }
    let eta: f64 = eta_numerator/eta_denominator;

    // FIND COUPLING MULTIPLIER
    // Split coupling multiplier into discrete chuncks [sic]->easier debugging
    let param1 = consts::CBET_CONST/(consts::OMEGA*(consts::TE_EV/1e3 + 3.0 * consts::TI_EV/1e3/consts::Z));
    let param2 = ne_over_nc/consts::IAW*consts::IAW*consts::IAW*eta; // Need to fix iaw
    let param3 = (eta*eta-1.0).powi(2) + consts::IAW*consts::IAW*eta*eta;
    let param4 = interaction_mult;
    // WILLOW SAYS: "ds" is replaced with crossing.dkmag, cuz that's what it is
    let coupling_mult = param1*param2/param3*param4*crossing.dkmag; // get the coupling multiplier
    
    let other_intensity1 = raycross.i_b;
    let other_intensity2 = raycross_next.i_b;
    // average the intensity of the other ray across the cell
    let avg_intensity = (other_intensity1+other_intensity2)/2.0;
    
    coupling_mult*avg_intensity
}

/// Initialize the i_b and w_mult values of each crossing, given an initial intensity. One TODO
/// would be to have the intensity be a property of the beam instead, so different beams can
/// have different intensities in the future.
///
/// C++ equivalent: initArrays in cbet.cpp
pub fn init_crossings(beams: &mut [Beam]) /*-> f64*/ {
    //let mut dkmags = Vec::new();
    beams.iter().for_each(|beam| {
        let increment = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(beam.rays.len() as f64-1.0);
        beam.rays.iter().enumerate().for_each(|(i, ray)| {
            let offset = consts::BEAM_MIN_Z + (increment*i as f64);
            let intensity = (beam.intensity/1e14)*f64::exp(-2.0*f64::abs(offset/2e-4).powf(4.0));

            ray.crossings.iter().for_each(|crossing| {
                let mut crossing = crossing.lock().unwrap();
                crossing.i_b = intensity;
                crossing.w_mult = 1.0;

                //if f64::abs(crossing.dkmag) > 1e-10 {
                //    dkmags.push(crossing.dkmag);
                //}
            });
        });
    });

    //dkmags.sort_unstable_by(|a, b| f64::partial_cmp(a, b).unwrap());
    //dkmags[dkmags.len()/2]
}

/// Populates the raystore vector of each beam, which stores a single ray number (as well as
/// a boolean for whether or not there is a ray there) per beam per spot on the mesh that the
/// beam crosses through.
fn create_raystore(beams: &mut [Beam], nx: usize, nz: usize) { 
    // beams, beams, they're good for your heart...
    beams.iter_mut().for_each(|beam| {
        beam.raystore = (0..nx*nz).map(|i| {
            let x = i / nz;
            let z = i % nz;
            let marked = &beam.marked[x*nz+z];
            match marked.len() {
                0 => (false, (0, 0)),
                1 => (true, marked[0]),
                // is it really worth it to use random for this?
                x => (true, marked[rand::random::<usize>() % x]),
            }
        }).collect();
    });
}
