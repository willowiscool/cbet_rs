use crate::consts;
use std::sync::Mutex;

/// TODO: Give beam a position and direction within it, and then make a general new beam
/// function that takes position and direction and populates rays. For now, the functions that
/// create the rays within the beams are hardcoded.
///
/// * rays: you can tell
/// * marked: a list of each spot on the mesh. each list contains the rays that pass through
/// it, stored by their index. It would perhaps be more correct to store the rays in as
/// `Rc<RefCell<Ray>>`, so that marked can own references to the specific rays, and they can
/// still be modified (as crossings are added), but I assume this would also make ray tracing
/// slower, so I think storing the indices is okay.
/// * raystore: similar to marked, except it only stores a single ray per each spot on the mesh
/// (by its index), as well as a boolean that notes whether or not there is a ray in that spot
/// at all. This is populated by cbet, which chooses a random ray when there are multiple rays
/// in the same zone.
#[derive(Debug)]
pub struct Beam {
    pub rays: Vec<Ray>,
    pub marked: Vec<Mutex<Vec<(usize, usize)>>>,
    pub raystore: Vec<(bool, (usize, usize))>,
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
                };
                ray.cz0 = ray.z0 + consts::CHILD_OFFSET;
                // c++ sets the first element of uray to interp(pow_x, phase_x + offset1, z0(n));
                // phase_x + offset1 is the same as z0
                // pow_x is a transformation of phase_x, so we can just
                // put in that calculation for pow_x
                // interp will just fetch the value of pow_x at z0
                // this might be wrong! so please correct it if it is!!
                // Also, I don't know why the pow stuff is like that haha
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
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
                };
                ray.cx0 = ray.x0 + consts::CHILD_OFFSET;
                // in this case, phase_x + offset2 = ray.x0
                ray
            }).collect(),
            marked: Vec::new(),
            raystore: Vec::new(),
        };
        beam
    }
}

/// Ray struct stores:
/// * List of crossings
/// * x0, z0: initial position
/// * cx0, cz0: initial position of child ray. the c++ impl. moves each ray in a single beam by
/// a child offset, so maybe it is better to have child offset x/z fields in the beam rather
/// than in the ray? But I'm defining it like this for now.
/// * kx0, kz0: initial velocity
/// * uray: absorbtion value or something. the c++ impl. uses a list for this but currently
/// only the first item of the list is stored so I'm only using one
#[derive(Debug)]
pub struct Ray {
    pub crossings: Vec<Crossing>,
    pub x0: f64,
    pub z0: f64,
    pub cx0: f64,
    pub cz0: f64,
    pub kx0: f64,
    pub kz0: f64,
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
#[derive(Debug)]
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

