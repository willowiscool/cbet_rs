pub mod consts;
pub mod mesh;
pub mod beam;
pub mod ray_trace;
pub mod utils;
pub mod cbet;
//pub use crate::rayTrace;

use crate::mesh::*;
use crate::beam::*;
use num::complex::Complex;

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

    // map returns field, has side effects into wplot
    let field = beams.iter().map(|beam| {
        let mut intensities = vec![0.0; mesh.nx*mesh.nz];
        let mut ct = vec![0.0; mesh.nx*mesh.nz];

        // create new field for each beam. code only uses the phase of ONE crossing in a zone
        // for the field reconstruction-- here we choose the last one.
        // an alternative could be to do the first one and skip the rest of the calculations but
        // i'm not sure if that's super necessary
        let mut field = vec![Complex::new(0.0,0.0); mesh.nx*mesh.nz];
        beam.rays.iter().for_each(|ray| {
            ray.crossings.iter().for_each(|crossing| {
                intensities[crossing.boxesx*mesh.nz + crossing.boxesz] += crossing.i_b;
                ct[crossing.boxesx*mesh.nz + crossing.boxesz] += 1.0;

                field[crossing.boxesx*mesh.nz + crossing.boxesz] = crossing.i_b * Complex::cis(crossing.phase);
            });
        });
        intensities.iter().enumerate().for_each(|(i, intensity)| {
            wplot[i] += (intensity / f64::max(ct[i], 1.0)).powi(2);
        });

        field
    }).reduce(|acc, x| {
        // sum fields
        // shout out https://www.reddit.com/r/rust/comments/eqc0hv/question_clean_way_to_do_elementwise_operations/
        std::iter::zip(acc, x).map(|(a, b)| a + b).collect::<Vec<Complex<f64>>>()
    }).unwrap();

    wplot = wplot.iter().map(|n| f64::sqrt(*n)).collect();
    let norm_field: Vec<f64> = field.iter().map(|n| n.norm()).collect();

    let wplot_dataset = file.new_dataset::<f64>()
        .shape([mesh.nx, mesh.nz])
        .create("wplot")?;
    wplot_dataset.write(
        &ndarray::Array::from_shape_vec((mesh.nx, mesh.nz), wplot)?
    )?;

    let field_dataset = file.new_dataset::<f64>()
        .shape([mesh.nx, mesh.nz])
        .create("field")?;
    field_dataset.write(
        &ndarray::Array::from_shape_vec((mesh.nx, mesh.nz), norm_field)?
    )?;

    Ok(())
}
