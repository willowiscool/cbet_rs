use std::sync::Mutex;
use ndarray::{Array3, Array2, Array};
use num::traits::Float;
use crate::consts;
use crate::mesh::XYZ;

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
    pub marked: Array3<Mutex<Vec<(usize, usize)>>>,
    pub raystore: Vec<(bool, (usize, usize))>,
    pub intensity: f64,
}
pub enum Source {
    MinX,
    MinY,
    MinZ,
}
impl Beam {
    /// Creates a new beam
    pub fn new(
        source: Source, offset: XYZ<f64>, nrays: usize, direction: XYZ<f64>, width: f64, intensity: f64
    ) -> Beam {
        // make rays: first, make circle of coordinates
        let coords = Array::linspace(-0.5, 0.5, nrays);
        Beam {
            rays: Array2::from_shape_fn((nrays, nrays), |(x, y)| (coords[x], coords[y]))
                .iter()
                .filter(|(x, y)| f64::sqrt(x.powi(2) + y.powi(2)) <= 0.5)
                .map(|(x, y)| {
                    let mut ray = Ray {
                        crossings: Vec::new(),
                        pos0: match source {
                            Source::MinX => XYZ {
                                x: consts::XMIN,
                                y: x * width + offset.y,
                                z: y * width + offset.z,
                            },
                            Source::MinY => XYZ {
                                x: x * width + offset.x,
                                y: consts::YMIN,
                                z: y * width + offset.z,
                            },
                            Source::MinZ => XYZ {
                                x: x * width + offset.x,
                                y: y * width + offset.y,
                                z: consts::ZMIN,
                            },
                        },
                        cpos0: XYZ { x: 0.0, y: 0.0, z: 0.0 },
                        k0: direction,
                    };
                    // GUESSING HERE!!!
                    ray.cpos0 = match source {
                        Source::MinX => XYZ {
                            x: consts::XMIN,
                            y: ray.pos0.y + consts::CHILD_OFFSET,
                            z: ray.pos0.z + consts::CHILD_OFFSET,
                        },
                        Source::MinY => XYZ {
                            x: ray.pos0.x + consts::CHILD_OFFSET,
                            y: consts::YMIN,
                            z: ray.pos0.z + consts::CHILD_OFFSET,
                        },
                        Source::MinZ => XYZ {
                            x: ray.pos0.x + consts::CHILD_OFFSET,
                            y: ray.pos0.y + consts::CHILD_OFFSET,
                            z: consts::ZMIN,
                        },
                    };
                    ray
                }).collect(),
            marked: Array3::from_shape_vec((0, 0, 0), vec![]).unwrap(),
            raystore: Vec::new(),
            intensity,
        }
    }
    /*
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
            marked: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
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
            marked: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
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
                let shift_amount = (ray.z0 - (consts::ZMIN + 1.0e-10)) / kz_norm;
                ray.x0 -= shift_amount * kx_norm;
                ray.z0 -= shift_amount * kz_norm;
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
            marked: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
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
                let shift_amount = (ray.z0 - (consts::ZMIN + 1.0e-10)) / kz_norm;
                ray.x0 -= shift_amount * kx_norm;
                ray.z0 -= shift_amount * kz_norm;
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
            marked: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
            raystore: Vec::new(),
            intensity: 5e15,
        };
        beam
    }*/
}

/// Ray struct stores:
/// * List of crossings (in mutexes for thread safety while referenced in Crossing!)
/// * pos0: initial position
/// * cpos0: initial position of child ray. the c++ impl. moves each ray in a single beam by
/// a child offset, so maybe it is better to have child offset x/z fields in the beam rather
/// than in the ray? But I'm defining it like this for now.
/// * k0: initial velocity
#[derive(Debug)]
pub struct Ray {
    pub crossings: Vec<Crossing>,
    pub pos0: XYZ<f64>,
    pub cpos0: XYZ<f64>,
    pub k0: XYZ<f64>,
}

/// Crossing struct stores
/// * pos: the real-world coordinates of the crossing (c++: crossesx/z)
/// * mesh_pos: the zone coordinates of the crossing (note: in c++, this is stored as
/// just boxes, where each element is a tuple)
/// * areaRatio: the ratio between parent and child ray
/// * dk, dkmag: vectors of movement to next crossing, i guess. computed in main.cpp after
/// ray tracing runs, but figured it would be more efficient to compute them in the ray tracing
/// function. also, computed as dkx_new, dkz_new, dkmag_new in cpp impl.
/// * i_b: the intensity, calculated in the CBET stage.
/// * energy: the energy multiplier from absorption calculations, multiplied into the initial
/// intensity --- exp(kds)!
/// * absorption_coeff: next energy/current energy
#[derive(Debug, Clone)]
pub struct Crossing {
    pub pos: XYZ<f64>,
    pub mesh_pos: XYZ<usize>,
    pub area_ratio: f64,
    pub dk: XYZ<f64>,
    pub dkmag: f64,
    pub i_b: f64,
    pub energy: f64,
    pub absorption_coeff: f64,
    pub phase: f64,
}
