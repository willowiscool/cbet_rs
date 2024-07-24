// use crate::consts;
use ndarray::Array2;
/// Point stores all of the information used at a given point in the mesh
///
/// kib_multiplier = 1e4 * kib * (eden / ncrit).powi(2) / f64::sqrt(dielectric)
/// where dielectric = 1 - eden / ncrit
///
/// permittivity_multiplier = permittivity * consts::OMEGA / consts::C_SPEED
#[derive(Debug)]
pub struct Point {
    pub x: f64,
    pub z: f64,
    pub eden: f64,
    pub machnum: f64,
    pub kib: f64,
    pub kib_multiplier: f64,
    pub permittivity_multiplier: f64,
}
/// Mesh struct
/// * dx, dz are width and height of each zone
/// * nx, nz are the number of zones in each dimension
/// * xmin, xmax, zmin, zmax store the minimum real-world coordinates (inclusive)
pub struct Mesh {
    pub points: Array2<Point>,
    pub dx: f64,
    pub dz: f64,
    pub nx: usize,
    pub nz: usize,
    pub xmin: f64,
    pub xmax: f64,
    pub zmin: f64,
    pub zmax: f64,
}

impl Mesh {
    /// Creates a new linear mesh. pointFn must take real x, z coordinates (not mesh coordinates)
    /// and return a Point.
    ///
    /// Consider changing the params to a Config struct?
    pub fn new_lin<F>(
        xmin: f64, xmax: f64, nx: usize, zmin: f64, zmax: f64, nz: usize,
        point_fn: F,
    ) -> Mesh
    where
        F: Fn(f64, f64) -> Point
    {
        let dx = (xmax - xmin)/(nx as f64 - 1.0);
        let dz = (zmax - zmin)/(nz as f64 - 1.0);
        Mesh {
            xmin,
            xmax,
            zmin,
            zmax,
            dx,
            dz,
            nx,
            nz,
            points: Array2::from_shape_fn(
                (nx, nz),
                |(meshx, meshz)| {
                    point_fn(meshx as f64 * dx + xmin, meshz as f64 * dz + zmin)
                }
            )
        }
    }

    /// Gets the point at mesh coordinates x, z (TODO: error check?)
    pub fn get(&self, x: usize, z: usize) -> &Point {
        //&self.points[(x*self.nz + z) as usize]
        self.points.get((x, z)).expect("Tried to get point out of bounds")
    }
    /// Gets the point at mesh coordinates x, z, mutably
    pub fn get_mut(&mut self, x: usize, z: usize) -> &mut Point {
        //&mut self.points[(x*self.nz + z) as usize]
        self.points.get_mut((x, z)).expect("Tried to get point out of bounds")
    }
    /// Gets the mesh coordinates of a point given its real-world coordinates.
    /// This code may be improved to be faster but it isn't used often so it's ok.
    /// I also think this assumes a linear mesh.
    /// Corresponds to lines 435-449 of launch_ray_XZ.cpp.
    /// TODO: Ideally, this would return a Result in case the coordinate given is out of
    /// bounds. Right now, it would just return (0.0, 0.0)
    pub fn get_mesh_coords(&self, x: f64, z: f64) -> (usize, usize) {
        self.get_mesh_coords_in_area(x, z, (0, 0), (self.nx-1, self.nz-1))
    }

    /// Gets the mesh coordinates of a point given its real-world coordinates, only
    /// searching within the subsection of the mesh given. The subsection is given by two
    /// (x, z) tuples - minpt and maxpt - and is inclusive. Corresponds to lines 498-511 of
    /// launch_ray_XZ.cpp
    pub fn get_mesh_coords_in_area(&self, x: f64, z: f64, minpt: (usize, usize), maxpt: (usize, usize)) -> (usize, usize) {
        let mut point = (0, 0);
        for xx in minpt.0..=maxpt.0 {
            let px = self.get(xx, 0).x;
            if (x - px <= (1.0 + 1.0e-10)*self.dx) &&
                (x - px >= -(0.0 + 1.0e-10)*self.dx) {
                point.0 = xx;
                break;
            }
        }
        for zz in minpt.1..=maxpt.1 {
            let pz = self.get(0, zz).z;
            if (z - pz <= (1.0 + 1.0e-10)*self.dz) &&
                (z - pz >= -(0.0 + 1.0e-10)*self.dz) {
                point.1 = zz;
                break;
            }
        }
        point
    }

    /* Commented: I don't think these values are used
    /// Returns the value, depending on eden, calculated at line 388 of main.cpp in the ray
    /// tracing function. Figured instead of calculating the values beforehand, just move
    /// it inline. Is this a bad move? TODO: perhaps, have wpe be a value of the mesh, so
    /// it can either be calculated for all eden values when the mesh is created or at
    /// least saved so multiple calculations of the same wpe value don't occur.
    pub fn wpe(&self, x: usize, z: usize) -> f64 {
        (self.get(x, z).eden*1e6*consts::EC*consts::EC/(consts::ME*consts::E0)).sqrt()
    }*/
}
