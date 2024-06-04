pub mod consts;
pub mod mesh;
pub mod beam;
pub mod ray_trace;
pub mod utils;
pub mod cbet;
//pub use crate::rayTrace;

use crate::mesh::*;
use crate::beam::*;

/// Creates a new mesh using the consts defined in consts.rs
pub fn new_default_mesh() -> Mesh {
    Mesh::new_lin(
        consts::XMIN, consts::XMAX, consts::NX,
        consts::ZMIN, consts::ZMAX, consts::NZ,
    )
}

/// Saves the data into an hdf5 file. Trying to be as identical to the c++ code as possible,
/// but in the future maybe saving the structs as data would be nicer?)
pub fn save_hdf5(mesh: &Mesh, beams: &[Beam], filename: &str) -> Result<(), hdf5::Error> {
    let file = hdf5::File::create(filename)?;
    // https://stackoverflow.com/questions/70568332/specifying-shape-of-dataset-when-creating-hdf5-file

    // named WPlot1 in c++
    let mut wplot = vec![0.0; mesh.nx*mesh.nz];
    beams.iter().for_each(|beam| {
        let mut intensities = vec![0.0; mesh.nx*mesh.nz];
        let mut ct = vec![0.0; mesh.nx*mesh.nz];
        beam.rays.iter().for_each(|ray| {
            ray.crossings.iter().for_each(|crossing| {
                intensities[crossing.boxesx*mesh.nz + crossing.boxesz] += crossing.i_b;
                ct[crossing.boxesx*mesh.nz + crossing.boxesz] += 1.0;
            });
        });
        intensities.iter().enumerate().for_each(|(i, intensity)| {
            wplot[i] += (intensity / f64::max(ct[i], 1.0)).powi(2);
        });
    });
    wplot = wplot.iter().map(|n| f64::sqrt(*n)).collect();

    let wplot_dataset = file.new_dataset::<f64>()
        .shape([mesh.nx, mesh.nz])
        .create("wplot")?;
    wplot_dataset.write(
        &ndarray::Array::from_shape_vec((mesh.nx, mesh.nz), wplot)?
    )?;

    Ok(())
}
