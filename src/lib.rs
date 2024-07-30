pub mod consts;
pub mod mesh;
pub mod beam;
pub mod ray_trace;
pub mod utils;
//pub mod cbet;
//pub use crate::rayTrace;

use crate::mesh::*;
use crate::beam::*;
use num::complex::Complex;
use ndarray::Array3;
use ndarray::Zip;

/// Creates a new mesh using the consts defined in consts.rs and here. It has a linear electron
/// density gradient. Based on main.cpp around line 500 or so
///
/// Not implemented as mesh::new since it's got way too many constants and stuff
pub fn new_mesh() -> Mesh {
    // constants
    let ncrit = consts::NCRIT;
    //let kib = 2.0 * consts::NU_EI_C * 1e12/(consts::C_SPEED*1e4);
    let kib = 0.0;
    let xmax = consts::XMAX;
    let xmin = consts::XMIN;
    Mesh::new_lin(
        XYZ { x: consts::XMIN, y: consts::YMIN, z: consts::ZMIN },
        XYZ { x: consts::XMAX, y: consts::YMAX, z: consts::ZMAX },
        XYZ { x: consts::NX, y: consts::NY, z: consts::NZ },
        |coord| {
            // sqrt_permittivity also = sqrt_dielectric
            //let eden = f64::max(0.0, ((0.4*ncrit-0.001*ncrit)/(xmax-xmin))*(x-xmin)+(0.001*ncrit));
            let eden = f64::max(0.0, ((0.4*ncrit-0.1*ncrit)/(xmax-xmin))*(coord.x-xmin)+(0.1*ncrit));
            let sqrt_permittivity = f64::sqrt(1.0 - (eden / ncrit));
            Point {
                pos: coord,
                eden,
                machnum: f64::min(0.0, (((-1.8)-(-1.0))/(xmax-xmin))*(coord.x-xmin))+(-1.0),
                kib,
                kib_multiplier: 1e4 * kib * (eden / ncrit).powi(2) / sqrt_permittivity,
                permittivity_multiplier: f64::max(sqrt_permittivity, 0.0) * consts::OMEGA / consts::C_SPEED,
            }
        }
    )
}

/// Saves the data into an hdf5 file. Trying to be as identical to the c++ code as possible,
/// but in the future maybe saving the structs as data would be nicer?)
pub fn save_hdf5(mesh: &Mesh, beams: &[Beam], filename: &str) -> Result<(), hdf5::Error> {
    let file = hdf5::File::create(filename)?;
    // https://stackoverflow.com/questions/70568332/specifying-shape-of-dataset-when-creating-hdf5-file

    // named WPlot1 in c++
    let mut wplot = Array3::zeros(mesh.n.tuple());

    // map returns field, has side effects into wplot
    let fields: Vec<Array3<Complex<f64>>> = beams.iter().map(|beam| {
        let mut intensities = Array3::zeros(mesh.n.tuple());
        let mut ct = Array3::zeros(mesh.n.tuple());

        // create new field for each beam. code only uses the phase of ONE crossing in a zone
        // for the field reconstruction-- here we choose the last one.
        // an alternative could be to do the first one and skip the rest of the calculations but
        // i'm not sure if that's super necessary
        let mut field = Array3::from_elem(mesh.n.tuple(), Complex::new(0.0, 0.0));
        beam.rays.iter().for_each(|ray| {
            ray.crossings.iter().for_each(|crossing| {
                *intensities.get_mut(crossing.mesh_pos.tuple()).unwrap() += crossing.i_b;
                *ct.get_mut(crossing.mesh_pos.tuple()).unwrap() += 1.0;

                *field.get_mut(crossing.mesh_pos.tuple()).unwrap() = crossing.i_b * Complex::cis(crossing.phase);
            });
        });
        // this is so snazzy
        Zip::from(&mut wplot)
            .and(&intensities)
            .and(&ct)
            .for_each(|wplot, &intensity: &f64, &ct| *wplot += (intensity / f64::max(ct, 1.0)).powi(2));

        field
    }).collect();
    let norm_field = Array3::from_shape_fn(
        mesh.n.tuple(),
        |(x, y, z)| fields.iter().fold(Complex::new(0.0, 0.0), |acc, field| acc + field[[x, y, z]]).norm()
    );

    wplot.map_inplace(|n| *n = f64::sqrt(*n));

    let wplot_dataset = file.new_dataset::<f64>()
        .shape([mesh.n.x, mesh.n.y, mesh.n.z])
        .create("wplot")?;
    wplot_dataset.write(&wplot)?;

    let field_dataset = file.new_dataset::<f64>()
        .shape([mesh.n.x, mesh.n.y, mesh.n.z])
        .create("field")?;
    field_dataset.write(&norm_field)?;

    Ok(())
}
