use std::sync::Mutex;
use crate::consts;

/// TODO: Give beam a position and direction within it, and then make a general new beam
/// function that takes position and direction and populates rays. For now, the functions that
/// create the rays within the beams are hardcoded.
///
/// * rays: you can tell
/// * marked: a list of each spot on the mesh. Each list contains the crossings that pass through
/// it, as an ray index and a crossing index.
/// * raystore: similar to marked, except it only stores a single crossing per each spot on the
/// mesh. It's an option because there may not be a crossing in that spot.
#[derive(Debug)]
pub struct Beam {
    pub rays: Vec<Ray>,
    pub marked: Vec<Vec<(usize, usize)>>,
    pub raystore: Vec<(bool, (usize, usize))>,
    pub intensity: f64,
}
impl Beam {
    /// Creates the first beam based on a bunch of constants in consts.rs. Close to 1:1
    /// copy of the c++ impl. in rayTracing.
    ///
    /// I choose to initialize the marked array in the ray tracing function, which gets the
    /// mesh (and thus the dimensions)
    pub fn beam1() -> Beam {
        let dz = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(consts::NRAYS as f64 - 1.0);
        let beam = Beam {
            rays: (0..consts::NRAYS).map(|i| {
                // can make uray with capacity NT but not sure that's needed
                let mut ray = Ray {
                    crossings: Vec::new(),
                    x0: consts::XMIN,
                    z0: i as f64 * dz + consts::BEAM_MIN_Z + consts::OFFSET1,
                    cx0: consts::XMIN,
                    cz0: 0.0,
                    kx0: consts::DIR1[0],
                    kz0: consts::DIR1[1],
                    uray: Vec::new(),
                };
                ray.cz0 = ray.z0 + consts::CHILD_OFFSET;
                // c++ sets the first element of uray to interp(pow_x, phase_x + offset1, z0(n));
                // phase_x + offset1 is the same as z0
                // pow_x is a transformation of phase_x, so we can just
                // put in that calculation for pow_x
                // interp will just fetch the value of pow_x at z0
                // this might be wrong! so please correct it if it is!!
                // Also, I don't know why the pow stuff is like that haha
                ray.uray.push(
                    consts::URAY_MULT*(-2.0*((ray.z0/consts::SIGMA).powi(2)).powf(4.0/2.0)).exp()
                );
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
            intensity: consts::INTENSITY,
        };
        beam
    }
    pub fn beam2() -> Beam {
        let dx = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(consts::NRAYS as f64 - 1.0);
        let beam = Beam {
            rays: (0..consts::NRAYS).map(|i| {
                // can make uray with capacity NT but not sure that's needed
                let mut ray = Ray {
                    crossings: Vec::new(),
                    x0: i as f64 * dx + consts::BEAM_MIN_Z + consts::OFFSET2,
                    z0: consts::ZMIN+0.1e-4,
                    cx0: 0.0,
                    cz0: consts::ZMIN+0.1e-4,
                    kx0: consts::DIR2[0],
                    kz0: consts::DIR2[1],
                    uray: Vec::new(),
                };
                ray.cx0 = ray.x0 + consts::CHILD_OFFSET;
                // in this case, phase_x + offset2 = ray.x0
                ray.uray.push(
                    consts::URAY_MULT2*(-2.0*((ray.x0/consts::SIGMA).powi(2)).powf(4.0/2.0)).exp()
                );
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
            intensity: consts::INTENSITY,
        };
        beam
    }
    /// This fn. is for the 2nd beam in the 3-beam test case (2d3beam branch of c++ impl.).
    /// The c++ impl. includes a "rotate" function that I'm attempting to recreate here
    pub fn beam3() -> Beam {
        let x = consts::OFFSET2;
        let z = consts::ZMIN + 0.1e-4;
        let kx = 0.6;
        let kz = 1.0;
        let k_mag = f64::sqrt(kx*kx + kz*kz);
        let kx_norm = kx / k_mag;
        let kz_norm = kz / k_mag;
        let dx = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(consts::NRAYS as f64 - 1.0);
        let dkx = (0.5286-0.4892)/(consts::NRAYS as f64 - 1.0);
        let dkz = (0.8670-0.8489)/(consts::NRAYS as f64 - 1.0);
        let beam = Beam {
            rays: (0..consts::NRAYS).map(|i| {
                let mut ray = Ray {
                    crossings: Vec::new(),
                    x0: (i as f64 * dx + consts::BEAM_MIN_Z) * (-kz_norm) + x,
                    z0: (i as f64 * dx + consts::BEAM_MIN_Z) * (kx_norm) + z,
                    cx0: 0.0,
                    cz0: 0.0,
                    kx0: i as f64 * dkx + 0.4892,
                    kz0: i as f64 * dkz + 0.8489,
                    uray: Vec::new(),
                };
                // ShiftZMin
                let shiftAmount = (ray.z0 - (consts::ZMIN + 1.0e-10)) / kz_norm;
                ray.x0 -= shiftAmount * kx_norm;
                ray.z0 -= shiftAmount * kz_norm;
                ray.uray.push(0.0);
                // rotateChild
                let curr_k_mag = f64::sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);
                assert!(curr_k_mag != 0.0);
                let curr_kx_norm = ray.kx0 / curr_k_mag;
                let curr_kz_norm = ray.kz0 / curr_k_mag;
                if curr_kx_norm > 0.0 {
                    ray.cx0 = consts::CHILD_OFFSET * (-curr_kz_norm) + ray.x0;
                    ray.cz0 = consts::CHILD_OFFSET * curr_kx_norm + ray.z0;
                } else {
                    ray.cx0 = consts::CHILD_OFFSET * curr_kz_norm + ray.x0;
                    ray.cz0 = consts::CHILD_OFFSET * (-curr_kx_norm) + ray.z0;
                }
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
            intensity: 5e15,
        };
        beam
    }
    /// This fn. is for the 3rd beam in the 3beam test case.
    pub fn beam4() -> Beam {
        let x = consts::OFFSET3;
        let z = consts::ZMIN;
        let kx = -0.3;
        let kz = 1.0;
        let k_mag = f64::sqrt(kx*kx + kz*kz);
        let kx_norm = kx / k_mag;
        let kz_norm = kz / k_mag;
        let dx = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/(consts::NRAYS as f64 - 1.0);
        let dkx = (-0.2971-(-0.2743))/(consts::NRAYS as f64 - 1.0);
        let dkz = (0.9548-0.9616)/(consts::NRAYS as f64 - 1.0);
        let beam = Beam {
            rays: (0..consts::NRAYS).map(|i| {
                let mut ray = Ray {
                    crossings: Vec::new(),
                    x0: (i as f64 * dx + consts::BEAM_MIN_Z) * (-kz_norm) + x,
                    z0: (i as f64 * dx + consts::BEAM_MIN_Z) * (kx_norm) + z,
                    cx0: 0.0,
                    cz0: 0.0,
                    kx0: i as f64 * dkx + -0.2743,
                    kz0: i as f64 * dkz + 0.9616,
                    uray: Vec::new(),
                };
                // ShiftZMin
                let shiftAmount = (ray.z0 - (consts::ZMIN + 1.0e-10)) / kz_norm;
                ray.x0 -= shiftAmount * kx_norm;
                ray.z0 -= shiftAmount * kz_norm;
                ray.uray.push(0.0);
                // rotateChild
                let curr_k_mag = f64::sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);
                assert!(curr_k_mag != 0.0);
                let curr_kx_norm = ray.kx0 / curr_k_mag;
                let curr_kz_norm = ray.kz0 / curr_k_mag;
                if curr_kx_norm > 0.0 {
                    ray.cx0 = consts::CHILD_OFFSET * (-curr_kz_norm) + ray.x0;
                    ray.cz0 = consts::CHILD_OFFSET * curr_kx_norm + ray.z0;
                } else {
                    ray.cx0 = consts::CHILD_OFFSET * curr_kz_norm + ray.x0;
                    ray.cz0 = consts::CHILD_OFFSET * (-curr_kx_norm) + ray.z0;
                }
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
            intensity: 5e15,
        };
        beam
    }
}

/// Ray struct stores:
/// * List of crossings (in mutexes for thread safety while referenced in Crossing!)
/// * x0, z0: initial position
/// * cx0, cz0: initial position of child ray. the c++ impl. moves each ray in a single beam by
/// a child offset, so maybe it is better to have child offset x/z fields in the beam rather
/// than in the ray? But I'm defining it like this for now.
/// * kx0, kz0: initial velocity
/// * uray: absorbtion value or something. the c++ impl. uses a list for this but currently
/// only the first item of the list is stored so I'm only using one
#[derive(Debug)]
pub struct Ray {
    pub crossings: Vec<Mutex<Crossing>>,
    pub x0: f64,
    pub z0: f64,
    pub cx0: f64,
    pub cz0: f64,
    pub kx0: f64,
    pub kz0: f64,
    pub uray: Vec<f64>,
}

/// Crossing struct stores
/// * x, z: the real-world coordinates of the crossing (c++: crossesx/z)
/// * boxesx, boxesz: the zone coordinates of the crossing (note: in c++, this is stored as
/// just boxes, where each element is a tuple)
/// * areaRatio: the ratio between parent and child ray
/// * dkx, dkz, dkmag: vectors of movement to next crossing, i guess. computed in main.cpp after
/// ray tracing runs, but figured it would be more efficient to compute them in the ray tracing
/// function. also, computed as dkx_new, dkz_new, dkmag_new in cpp impl.
/// * i_b: the intensity, calculated in the CBET stage.
/// * wMult: "cumulative product of the normalized ray energies" - cbet.cpp:336
#[derive(Debug, Clone)]
pub struct Crossing {
    pub x: f64,
    pub z: f64,
    pub boxesx: usize,
    pub boxesz: usize,
    pub area_ratio: f64,
    pub dkx: f64,
    pub dkz: f64,
    pub dkmag: f64,
    pub i_b: f64,
    pub w_mult: f64,
}
