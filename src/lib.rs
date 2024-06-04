pub mod consts;
pub mod mesh;
pub mod beam;
pub mod ray_trace;
pub mod utils;
pub mod cbet;
//pub use crate::rayTrace;

use crate::mesh::*;

/// Creates a new mesh using the consts defined in consts.rs
pub fn new_default_mesh() -> Mesh {
    Mesh::new_lin(
        consts::XMIN, consts::XMAX, consts::NX,
        consts::ZMIN, consts::ZMAX, consts::NZ,
    )
}
