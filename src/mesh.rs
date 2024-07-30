use crate::consts;
use ndarray::Array3;
#[derive(Debug, Clone, Copy)]
pub struct XYZ<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}
impl<T: Copy> XYZ<T> {
    pub fn tuple(&self) -> (T, T, T) {
        (self.x, self.y, self.z)
    }
}
/// Point stores all of the information used at a given point in the mesh
///
/// kib_multiplier = 1e4 * kib * (eden / ncrit).powi(2) / f64::sqrt(dielectric)
/// where dielectric = 1 - eden / ncrit
///
/// permittivity_multiplier = permittivity * consts::OMEGA / consts::C_SPEED
#[derive(Debug)]
pub struct Point {
    pub pos: XYZ<f64>,
    pub eden: f64,
    pub machnum: f64,
    pub kib: f64,
    pub kib_multiplier: f64,
    pub permittivity_multiplier: f64,
}
/// Mesh struct
/// * d are width and height of each zone
/// * n are the number of zones in each dimension
/// * min, max store the minimum real-world coordinates (inclusive)
pub struct Mesh {
    pub points: Array3<Point>,
    pub d: XYZ<f64>,
    pub n: XYZ<usize>,
    pub min: XYZ<f64>,
    pub max: XYZ<f64>,
}

impl Mesh {
    /// Creates a new linear mesh. pointFn must take real x, z coordinates (not mesh coordinates)
    /// and return a Point.
    ///
    /// Consider changing the params to a Config struct?
    pub fn new_lin<F>(
        min: XYZ<f64>, max: XYZ<f64>, n: XYZ<usize>,
        point_fn: F,
    ) -> Mesh
    where
        F: Fn(XYZ<f64>) -> Point
    {
        let d = XYZ {
            x: (max.x - min.x)/(n.x as f64 - 1.0),
            y: (max.y - min.y)/(n.y as f64 - 1.0),
            z: (max.z - min.z)/(n.z as f64 - 1.0),
        };
        Mesh {
            min,
            max,
            d,
            n,
            points: Array3::from_shape_fn(
                n.tuple(),
                |(x, y, z)| {
                    point_fn(XYZ {
                        x: x as f64 * d.x + min.x,
                        y: y as f64 * d.y + min.y,
                        z: z as f64 * d.z + min.z,
                    })
                }
            )
        }
    }

    /// Gets the point at mesh coordinates x, z (TODO: error check?)
    pub fn get(&self, x: usize, y: usize, z: usize) -> &Point {
        //&self.points[(x*self.nz + z) as usize]
        self.points.get((x, y ,z)).expect("Tried to get point out of bounds")
    }
    /// Gets the point at mesh coordinates x, z, mutably
    pub fn get_mut(&mut self, x: usize, y: usize, z: usize) -> &mut Point {
        //&mut self.points[(x*self.nz + z) as usize]
        self.points.get_mut((x, y, z)).expect("Tried to get point out of bounds")
    }
    /// Gets the mesh coordinates of a point given its real-world coordinates.
    /// This code may be improved to be faster but it isn't used often so it's ok.
    /// I also think this assumes a linear mesh.
    /// Corresponds to lines 435-449 of launch_ray_XZ.cpp.
    /// TODO: Ideally, this would return a Result in case the coordinate given is out of
    /// bounds. Right now, it would just return (0.0, 0.0)
    pub fn get_mesh_coords(&self, x: f64, y: f64, z: f64) -> XYZ<usize> {
        self.get_mesh_coords_in_area(
            XYZ { x, y, z },
            XYZ { x: 0, y: 0, z: 0 },
            XYZ { x: self.n.x-1, y: self.n.y-1, z: self.n.z-1 }
        )
    }

    /// Gets the mesh coordinates of a point given its real-world coordinates, only
    /// searching within the subsection of the mesh given. The subsection is given by two
    /// (x, z) tuples - minpt and maxpt - and is inclusive. Corresponds to lines 498-511 of
    /// launch_ray_XZ.cpp
    pub fn get_mesh_coords_in_area(&self, coord: XYZ<f64>, minpt: XYZ<usize>, maxpt: XYZ<usize>) -> XYZ<usize> {
        let mut point = XYZ { x: 0, y: 0, z: 0 };
        for xx in minpt.x..=maxpt.x {
            let px = self.get(xx, 0, 0).pos.x;
            if (coord.x - px <= (1.0 + 1.0e-10)*self.d.x) &&
                (coord.x - px >= -(0.0 + 1.0e-10)*self.d.x) {
                point.x = xx;
                break;
            }
        }
        for yy in minpt.y..=maxpt.y {
            let py = self.get(0, yy, 0).pos.y;
            if (coord.y - py <= (1.0 + 1.0e-10)*self.d.y) &&
                (coord.y - py >= -(0.0 + 1.0e-10)*self.d.y) {
                point.y = yy;
                break;
            }
        }
        for zz in minpt.z..=maxpt.z {
            let pz = self.get(0, 0, zz).pos.z;
            if (coord.z - pz <= (1.0 + 1.0e-10)*self.d.z) &&
                (coord.z - pz >= -(0.0 + 1.0e-10)*self.d.z) {
                point.z = zz;
                break;
            }
        }
        point
    }

    /// Returns the value, depending on eden, calculated at line 388 of main.cpp in the ray
    /// tracing function. Figured instead of calculating the values beforehand, just move
    /// it inline. Is this a bad move? TODO: perhaps, have wpe be a value of the mesh, so
    /// it can either be calculated for all eden values when the mesh is created or at
    /// least saved so multiple calculations of the same wpe value don't occur.
    pub fn wpe_squared(&self, x: usize, y: usize, z: usize) -> f64 {
        self.get(x, y, z).eden*1e6*consts::EC*consts::EC/(consts::ME*consts::E0)
    }
}
