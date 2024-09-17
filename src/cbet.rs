use crate::beam::*;
use crate::mesh::*;
use crate::consts;

#[cxx::bridge]
mod ffi {
    // shared structs: fields visible to both languages
    #[derive(Clone)]
    pub struct CbetCrosses {
        pub b_num: usize,
        pub r_num: usize,
        pub c_num: usize,
        pub c_num_next: usize,
        pub coupling_mult: f64,
    }
    /// Minimal crossing struct that contains only the info needed for cbet
    #[derive(Clone)]
    pub struct CbetCrossing {
        pub intensity: f64,
        pub absorption_coeff: f64,
    }

    unsafe extern "C++" {
        include!("cbet_rs/include/cpp_cbet.cuh");

        unsafe fn cpp_cbet(cbet_crossings: *mut CbetCrossing, cbet_crosses: *mut CbetCrosses, nb: usize, nr: usize, nc: usize);
    }
}
use ffi::{CbetCrosses, CbetCrossing};

/// Do the CBET calculation! Populates the i_b fields of each crossing
pub fn cbet(mesh: &Mesh, beams: &mut [Beam]) {
    consts::init_consts(); // inits CS const
    create_raystore(beams, mesh.nx, mesh.nz);

    let (
        mut cbet_crossings,
        mut cbet_crosses,
        (nb, nr, nc),
    ) = create_cbet_crossings(mesh, beams);

    unsafe {
        ffi::cpp_cbet(cbet_crossings.as_mut_ptr(), cbet_crosses.as_mut_ptr(), nb, nr, nc);
    }
    println!("Finished cbet, abt to post");
    post(mesh, beams, &cbet_crossings, nr, nc);
}

/// Turn beam crossings into cbet_crossings which only contain the necessary info needed.
pub fn create_cbet_crossings(mesh: &Mesh, beams: &[Beam]) -> (Vec<CbetCrossing>, Vec<CbetCrosses>, (usize, usize, usize)) {
    // dimensions
    let nb = beams.len();
    let nr = consts::NRAYS;
    let nc = beams.iter().fold(0, |nc, beam| usize::max(nc,
        beam.rays.iter().fold(0, |nc, ray| usize::max(nc, ray.crossings.len()))
    ));

    let mut cbet_crossings = vec![CbetCrossing {
        intensity: 0.0,
        absorption_coeff: 0.0,
    }; nb*nr*nc];
    // this may be taking up too much memory in the future,
    // as in reality the vast majority of these structs (unlike the ones above)
    // will be empty, but this is a good first try.
    //
    // The goal of doing things this way is to have zero indirections,
    // so copying to GPU memory is easier.
    let mut cbet_crosses = vec![CbetCrosses {
        b_num: 0,
        r_num: 0,
        c_num: 0,
        c_num_next: 0,
        coupling_mult: 0.0,
    }; nb*nr*nc*nb];

    beams.iter().enumerate().for_each(|(bn, beam)| {
        beam.rays.iter().enumerate().for_each(|(rn, ray)| {
            ray.crossings.iter().enumerate().for_each(|(cn, crossing)| {
                cbet_crossings[
                    (((bn*nr)+rn)*nc)+cn
                ] = CbetCrossing {
                    intensity: crossing.i_b,
                    absorption_coeff: crossing.absorption_coeff,
                };
                let mut crosses_offset = 0;
                beams.iter().enumerate()
                    // o is for other
                    .filter(|(o_b_num, other_beam)|
                        *o_b_num != bn &&
                        other_beam.raystore[crossing.boxesx*mesh.nz+crossing.boxesz].0)
                    .for_each(|(o_b_num, o_beam)| {
                        // get raycross, raycross_next
                        let (_, (o_rayind, raycross_ind)) = o_beam.raystore[crossing.boxesx*mesh.nz+crossing.boxesz];
                        let o_ray_crossings = &o_beam.rays[o_rayind].crossings;
                        let raycross = &o_ray_crossings[raycross_ind];
                        let raycross_next_ind = usize::min(raycross_ind+1, o_ray_crossings.len()-1);
                        let raycross_next = &o_ray_crossings[raycross_next_ind];

                        cbet_crosses[
                            (((((bn*nr)+rn)*nc)+cn)*nb)+crosses_offset
                        ] = CbetCrosses {
                            b_num: o_b_num,
                            r_num: o_rayind,
                            c_num: raycross_ind,
                            c_num_next: raycross_next_ind,
                            coupling_mult: get_coupling_mult(mesh, crossing, raycross, raycross_next)
                        };
                        crosses_offset += 1;
                    })
            });
        });
    });
    (cbet_crossings, cbet_crosses, (nb, nr, nc))
}

/// Updates the i_b values of each crossing, one last time, while also moving this info from the
/// cbet_crossings back into the beams.
fn post(mesh: &Mesh, beams: &mut [Beam], cbet_crossings: &Vec<CbetCrossing>, nr: usize, nc: usize) {
    let w0 = 2.0*std::f64::consts::PI*consts::C_SPEED/consts::LAMBDA;
    let norm_factor_const = f64::sqrt(8.0*std::f64::consts::PI/consts::C_SPEED) * consts::ESTAT / (consts::ME_G*consts::C_SPEED*w0) * f64::sqrt(1e14*1e7);

    beams.iter_mut().enumerate().for_each(|(b_num, beam)| {
        beam.rays.iter_mut().enumerate().for_each(|(r_num, ray)| {
            // doing this instead of an iterator because there is an immutable
            // borrow of TWO crossings (to calculate area_avg) followed by a
            // mutable borrow of ONE crossing
            for c_num in 0..ray.crossings.len() {
                let crossing = &ray.crossings[c_num];
                let crossing_next = &ray.crossings[usize::min(c_num+1, ray.crossings.len()-1)];
                // some of this code looks copied from the cbet math code (TODO)
                let ix = crossing.boxesx;
                let iz = crossing.boxesz;

                let area_avg = (crossing.area_ratio+crossing_next.area_ratio)/2.0;
                let ne_over_nc = mesh.get(ix, iz).eden/consts::NCRIT;
                let ne_over_nc_corrected = f64::min(ne_over_nc, 1.0);
                let ne_term = f64::sqrt(1.0-ne_over_nc_corrected);
                let epsilon_eff = ne_term*ne_term;
                let interaction_mult = 1.0/(area_avg*ne_term)*1.0/f64::sqrt(epsilon_eff);
                let norm_factor = norm_factor_const * f64::sqrt(interaction_mult) * epsilon_eff.powf(0.25);
                let intensity_new = f64::sqrt(
                    cbet_crossings[(((b_num*nr)+r_num)*nc)+c_num].intensity
                ) * norm_factor;

                // now we have mutable borrow
                ray.crossings[c_num].i_b = intensity_new;
            }
        });
    });
}

/// This is where all of the math lives. Lines 162 through 209 of cbet.cpp are translated here.
/// Only thing that is changed is variable names are changed to snake case cuz I don't want
/// cargo to yell at me :(
fn get_coupling_mult(mesh: &Mesh, crossing: &Crossing, raycross: &Crossing, raycross_next: &Crossing) -> f64 {
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
    // random polarization (TODO include a flag for this?) taken from commit b44a of c++ code
    // coupling_mult *= (1.0 + (kx_seed * kx_pump + kz_seed * kz_pump).powi(2)) / 4.0;
    coupling_mult
}

/// Initialize the i_b and w_mult values of each crossing, given an initial intensity.
///
/// C++ equivalent: initArrays in cbet.cpp
pub fn init_crossings(beams: &mut [Beam]) /*-> f64*/ {
    //let mut dkmags = Vec::new();
    beams.iter_mut().for_each(|beam| {
        let increment = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(beam.rays.len() as f64-1.0);
        beam.rays.iter_mut().enumerate().for_each(|(i, ray)| {
            let offset = consts::BEAM_MIN_Z + (increment*i as f64);
            let intensity = (beam.intensity/1e14)*f64::exp(-2.0*f64::abs(offset/2e-4).powf(4.0));

            for j in 0..ray.crossings.len() {
                ray.crossings[j].i_b = intensity * ray.crossings[j].energy;
                /*if j == ray.crossings.len()-1 {
                    ray.crossings[j].absorption_coeff = 1.0;
                } else {
                    ray.crossings[j].absorption_coeff = ray.crossings[j+1].energy / ray.crossings[j].energy;
                }*/
                if j == 0 {
                    ray.crossings[j].absorption_coeff = ray.crossings[j].energy;
                } else {
                    ray.crossings[j].absorption_coeff = ray.crossings[j].energy / ray.crossings[j-1].energy;
                }
            }
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
            let marked = &beam.marked[[x, z]].lock().unwrap();
            match marked.len() {
                0 => (false, (0, 0)),
                1 => (true, marked[0]),
                // is it really worth it to use random for this?
                x => (true, marked[rand::random::<usize>() % x]),
            }
        }).collect();
    });
}
