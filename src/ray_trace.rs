use crate::beam::*;
use crate::mesh::*;
use crate::utils;
use crate::consts;
use std::cmp::{min, max};
use std::sync::Mutex;
use rayon::prelude::*;
use ndarray::Array3;

/// The rayTrace function populates the crossings of the rays in the beam array. It also
/// populates the marked array of each beam. It overwrites any existing crossings and marked
/// arrays.
///
/// It relies heavily on constants defined in consts.rs, under `// Ray tracing constants`.
/// TODO: stop doing that.
pub fn ray_trace(mesh: &Mesh, beams: &mut [Beam]) {
    // caclulate dedendz, dedendx
    // dedendx/z are being put into one vec which holds an (x, z) tuple
    // deden = derivative of eden
    // need to skip points on outer edges of mesh

    let mut deden = Array3::from_shape_fn(mesh.points.dim(), |(x, y, z)| {
        if x == mesh.n.x - 1 || y == mesh.n.y - 1 || z == mesh.n.z - 1 {
            XYZ { x: 0.0, y: 0.0, z: 0.0 }
        } else {
            let point = mesh.get(x, y, z);
            XYZ {
                x: if x == 0 { 0.0 } else {
                    let prev_point_x = mesh.get(x-1, y, z);
                    (point.eden - prev_point_x.eden)/(point.pos.x-prev_point_x.pos.x)
                },
                y: if y == 0 { 0.0 } else {
                    let prev_point_y = mesh.get(x, y-1, z);
                    (point.eden - prev_point_y.eden)/(point.pos.y-prev_point_y.pos.y)
                },
                z: if z == 0 { 0.0 } else {
                    let prev_point_z = mesh.get(x, y, z-1);
                    (point.eden - prev_point_z.eden)/(point.pos.z-prev_point_z.pos.z)
                },
            }
        }
    });

    // probably can optimize loop order but this is only done once so it's fine
    for i in 0..mesh.n.x {
        for j in 0..mesh.n.y {
            // for every z = 0, z = mesh.n.z - 1
            deden[[i, j, 0]].z = deden[[i, j, 1]].z;
            // XYZ implements copy so no need to worry
            deden[[i, j, mesh.n.z-1]] = deden[[i, j, mesh.n.z-2]];
        }
        for j in 0..mesh.n.z {
            // for every y = 0, y = mesh.n.y - 1
            deden[[i, 0, j]].y = deden[[i, 1, j]].y;
            deden[[i, mesh.n.y-1, j]] = deden[[i, mesh.n.y-2, j]];
        }
    }
    for i in 0..mesh.n.y {
        for j in 0..mesh.n.z {
            deden[[0, i, j]].x = deden[[1, i, j]].x;
            deden[[mesh.n.x-1, i, j]] = deden[[mesh.n.x-2, i, j]];
        }
    }

    beams.iter_mut().for_each(|beam| {
        // Initialize the marked array
        //beam.marked = Array2::from_shape_fn((mesh.nx, mesh.nz), |_| Mutex::new(Vec::new()));
        beam.marked = Array3::from_shape_fn((mesh.n.x, mesh.n.y, mesh.n.z), |_| Mutex::new(Vec::new()));
        beam.rays.par_iter_mut().enumerate().for_each(|(i, ray)| {
            let child = launch_child_ray(ray, mesh, &deden);
            launch_parent_ray(ray, mesh, &deden, &child, &beam.marked, i);
        });
    });
}

/// Identical code for parent/child so I moved it out. Gets the k constant that is used to
/// calculate the initial velocities
fn get_k(mesh: &Mesh, mesh_coords: XYZ<usize>) -> f64 {
    // NOTE: it looks like this assumes a linear electron density gradient!
    // The Mesh::wpe function (currently commented out) finds the wpe without doing that,
    // but only for the point, not interpolated. Consider using that to find wpe instead?
    /*let wpe_interp = f64::sqrt(
        utils::interp_closure(
            |i| { mesh.points[[i, 0]].eden }, mesh.nx, // eden.col(0) as y, ysize
            |i| { mesh.points[[i, 0]].x }, mesh.nx, // x.col(0) as x, xsize
            x0
        ) * 1e6 * consts::EC*consts::EC / (consts::ME*consts::E0)
    );*/
    let wpe_squared = mesh.wpe_squared(mesh_coords.x, mesh_coords.y, mesh_coords.z);
    // this is k
    f64::sqrt((consts::OMEGA.powi(2) - wpe_squared) / (consts::C_SPEED.powi(2)))
}

/// Launch the parent ray. Updates the marked vector, as well as the ray's crossings vector,
/// which it overwrites. More detailed comments for some of the code may be in
/// `launch_child_ray`, which I wrote first.
///
/// i is the ray index, because rays are stored in the marked array as indices (rather than
/// pointers, which might be what we could switch it to?)
fn launch_parent_ray(
    ray: &mut Ray, mesh: &Mesh, deden: &Array3<XYZ<f64>>,
    (childx, childy, childz): &(Vec<f64>, Vec<f64>, Vec<f64>), marked: &Array3<Mutex<Vec<(usize, usize)>>>, raynum: usize
) {
    // Maybe it is better to have all of these in one vector that stores a struct, i.e.
    // struct Timestamp { x, z, vx, vz, mx, mz }??
    // but I'll keep it like this for now, why not

    // looks like we only need these: + vx, vz later on
    let mut x = ray.pos0.x;
    let mut y = ray.pos0.y;
    let mut z = ray.pos0.z;

    // convert child coordinates into a set of distances (cumulative)
    // TODO: possibly move this code to launch_child_ray and do it in place?
    let mut distance = Vec::with_capacity(childx.len());
    distance.push(0.0);
    for i in 1..childx.len() {
        distance.push(distance[i-1] + f64::sqrt(
            (childx[i] - childx[i-1]).powi(2) +
            (childy[i] - childy[i-1]).powi(2) +
            (childz[i] - childz[i-1]).powi(2)
        ));
    }

    let init_area = {
        let init_diff_x = ray.pos0.x - childx[0];
        let init_diff_y = ray.pos0.y - childy[0];
        let init_diff_z = ray.pos0.z - childz[0];
        let init_diff_mag = f64::sqrt(init_diff_x.powi(2) + init_diff_y.powi(2) + init_diff_z.powi(2));
        let init_proj_coeff = f64::abs(ray.kx0 * (init_diff_z/init_diff_mag) - ray.kz0 * (init_diff_x/init_diff_mag));
        init_diff_mag*init_proj_coeff
    };

    // renamed from thisx_0, thisz_0, as well as meshx, meshz
    let (mut meshx, mut meshz) = mesh.get_mesh_coords(ray.x0, ray.z0);

    let k = get_k(mesh, ray.x0);
    let knorm = f64::sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);

    let mut vx = consts::C_SPEED*consts::C_SPEED * ((ray.kx0 / knorm) * k) / consts::OMEGA;
    let mut vz = consts::C_SPEED*consts::C_SPEED * ((ray.kz0 / knorm) * k) / consts::OMEGA;

    let mut kds = 0.0;
    let mut phase = 0.0;

    let mut curr_dist = 0.0;
    for _ in 1..consts::NT {
        // 1. Calculate velocity and position at current timestamp
        // ===
        // **copied from launch_child_ray**
        let dedendx = if meshx == mesh.nx - 1 { deden[[meshx, meshz]].0 } else {
            deden[[meshx, meshz]].0 + (
                (deden[[meshx+1, meshz]].0 - deden[[meshx, meshz]].0)
                / (mesh.get(meshx+1, meshz).x - mesh.get(meshx, meshz).x)
                * (x - mesh.get(meshx, meshz).x)
            )
        };

        // update pos, vel, equivalent to setting myx/z[tt] and myvx/z[tt]
        let prev_vx = vx;
        let prev_vz = vz;
        vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * dedendx * consts::DT;
        vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * deden[[meshx, meshz]].1 * consts::DT;
        let prev_x = x;
        let prev_z = z;
        x = prev_x + vx * consts::DT;
        z = prev_z + vz * consts::DT;

        // 2. update meshx, meshz
        let prev_meshx = meshx;
        let prev_meshz = meshz;
        (meshx, meshz) = mesh.get_mesh_coords_in_area(
            x, z,
            ( // min pt
                if meshx <= 1 { 0 } else { meshx - 1 },
                if meshz <= 1 { 0 } else { meshz - 1 },
            ),
            (min(mesh.nx - 1, meshx + 1), min(mesh.nz - 1, meshz + 1))
        );
        let meshpt = mesh.get(meshx, meshz);

        // defining an inline function since it takes so many variables from its
        // environment...
        let new_crossing = |x: f64, z: f64, frac: f64| {
            let distance_to_crossing = f64::sqrt((x - prev_x).powi(2) + (z - prev_z).powi(2));
            Crossing {
                x,
                z,
                boxesx: meshx,
                boxesz: meshz,
                area_ratio: {
                    let childxp = utils::interp(childx, &distance, curr_dist + distance_to_crossing);
                    let childzp = utils::interp(childz, &distance, curr_dist + distance_to_crossing);

                    let diff_x = x - childxp;
                    let diff_z = z - childzp;
                    let diff_mag = f64::sqrt(diff_x*diff_x + diff_z*diff_z);

                    let interpkx = frac*prev_vx + (1.0-frac)*vx;
                    let interpkz = frac*prev_vz + (1.0-frac)*vz;
                    let interpk_mag = f64::sqrt(interpkx*interpkx + interpkz*interpkz);

                    let proj_coeff = f64::abs(
                        (interpkx/interpk_mag) * (diff_z/diff_mag)
                        - (interpkz/interpk_mag) * (diff_x/diff_mag));

                    diff_mag*proj_coeff/init_area
                },
                dkx: 0.0,
                dkz: 0.0,
                dkmag: 0.0,
                i_b: -1.0,
                energy: f64::exp(kds - distance_to_crossing * meshpt.kib_multiplier),
                absorption_coeff: 0.0,
                phase: phase + distance_to_crossing * meshpt.permittivity_multiplier,
            }
        };

        let mut lastx = 10000.0;
        let mut lastz = 10000.0;
        let mut is_cross_x = false;
        let mut is_cross_z = false;

        // If crosses x, create a crossing
        // I changed the code quite a bit from the c++ impl, to eliminate redundancies
        // and make it more efficient and have less code... let's hope it actually
        // works haha
        // Looks like this assumes a linear mesh
        // NOTE: the modifications I've made assume no more than one crossing per
        // dimension in any given trajectory!!!
        if meshx != prev_meshx {
            // get real-world x coord of crossing
            // this MIGHT be mesh.get(max(meshx, prev_meshx), 0)
            // so TODO change that once you figure out this whole thing works
            let mut currx = mesh.get(meshx, 0).x;
            if !((x > currx && prev_x <= currx) || (x < currx && prev_x >= currx)) {
                currx = mesh.get(prev_meshx, 0).x;
            }
            assert!((x > currx && prev_x <= currx) || (x < currx && prev_x >= currx));

            // find the z point of intersection
            // TODO inline
            // don't let the name fool you, crossx is a z coordinate
            let crossx = utils::interp(&vec![prev_z, z], &vec![prev_x, x], currx);
            let frac = (currx - prev_x) / (x - prev_x);
            assert!((frac >= 0.0) && (frac <= 1.0));

            if f64::abs(crossx - lastz) > 1.0e-12 {
                let crossing = new_crossing(currx, crossx, frac);
                let mut markedpt = marked[[meshx, meshz]].lock().unwrap();
                markedpt.push((raynum, ray.crossings.len()));
                ray.crossings.push(crossing);

                is_cross_x = true;
                lastx = currx;
            }
        }

        // find z crossing :)
        if meshz != prev_meshz {
            // get real-world x coord of crossing
            // this MIGHT be mesh.get(max(meshx, prev_meshx), 0)
            // so TODO change that once you figure out this whole thing works
            let mut currz = mesh.get(0, meshz).z;
            if !((z > currz && prev_z <= currz) || (z < currz && prev_z >= currz)) {
                currz = mesh.get(0, prev_meshz).z;
            }
            assert!((z > currz && prev_z <= currz) || (z < currz && prev_z >= currz));

            // find the x point of intersection
            // TODO inline
            // don't let the name fool you, crossz is a x coordinate
            let crossz = utils::interp(&vec![prev_x, x], &vec![prev_z, z], currz);
            let frac = (currz - prev_z) / (z - prev_z);
            assert!((frac >= 0.0) && (frac <= 1.0));

            if f64::abs(crossz - lastx) > 1.0e-12 {
                let crossing = new_crossing(crossz, currz, frac);
                let mut markedpt = marked[[meshx, meshz]].lock().unwrap();
                markedpt.push((raynum, ray.crossings.len()));
                ray.crossings.push(crossing);

                is_cross_z = true;
                lastz = currz;
            }
        }
        // swap if out of order. then, calculate dkx, dkz, dkmag for prev. crossing(s)
        if is_cross_x && is_cross_z {
            let last_crossing_ind = ray.crossings.len()-1;
            let need_to_swap = {
                let last_crossing = &ray.crossings[last_crossing_ind];
                let second_last_crossing = &ray.crossings[last_crossing_ind-1];
                (x - prev_x) * (last_crossing.x - second_last_crossing.x) < 0.0
            };
            if need_to_swap {
                ray.crossings.swap(last_crossing_ind, last_crossing_ind-1);
            }
            if last_crossing_ind > 1 {
                calc_dk(ray, last_crossing_ind-2);
            }
            calc_dk(ray, last_crossing_ind-1);
        } else if (is_cross_x || is_cross_z) && ray.crossings.len() > 1 {
            calc_dk(ray, ray.crossings.len()-2);
        }
        let distance_travelled = f64::sqrt((x - prev_x).powi(2) + (z - prev_z).powi(2));
        curr_dist += distance_travelled;
        kds -= distance_travelled * meshpt.kib_multiplier;
        phase += distance_travelled * meshpt.permittivity_multiplier;
        if x < mesh.xmin || x > mesh.xmax || z < mesh.zmin || z > mesh.zmax {
            break;
        }
    }
}

/// Calculates and inserts dk values of crossing at ind. Ind cannot be the last index of the
/// crossings.
fn calc_dk(ray: &mut Ray, ind: usize) {
    assert!(ind < ray.crossings.len()-1);
    let dkx = ray.crossings[ind+1].x - ray.crossings[ind].x;
    let dkz = ray.crossings[ind+1].z - ray.crossings[ind].z;
    let dkmag = f64::sqrt(dkx.powi(2) + dkz.powi(2));
    // normalize, as in cpp impl.
    let dkx_new = dkx/dkmag;
    let dkz_new = dkz/dkmag;
    let dkmag_new = dkmag*10000.0;

    ray.crossings[ind].dkx = dkx_new;
    ray.crossings[ind].dkz = dkz_new;
    ray.crossings[ind].dkmag = dkmag_new;
}

/// Launch the child ray. Returns a vector of coordinates at each timestamp that represent the
/// trajectory of the child ray.
fn launch_child_ray(
    ray: &mut Ray, mesh: &Mesh, deden: &Array3<XYZ<f64>>
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();
    let mut z: Vec<f64> = Vec::new();
    x.push(ray.cpos0.x);
    y.push(ray.cpos0.y);
    z.push(ray.cpos0.z);

    // renamed from thisx_0, thisz_0
    let XYZ { x: mut meshx, y: mut meshy, z: mut meshz } = mesh.get_mesh_coords(ray.cpos0.x, ray.cpos0.y, ray.cpos0.z);

    let k = get_k(mesh, XYZ { x: meshx, y: meshy, z: meshz });
    let knorm = f64::sqrt(ray.k0.x.powi(2) + ray.k0.y.powi(2) + ray.k0.z.powi(2));
    let mut vx = consts::C_SPEED.powi(2) * ((ray.k0.x / knorm) * k) / consts::OMEGA;
    let mut vy = consts::C_SPEED.powi(2) * ((ray.k0.y / knorm) * k) / consts::OMEGA;
    let mut vz = consts::C_SPEED.powi(2) * ((ray.k0.z / knorm) * k) / consts::OMEGA;

    for tt in 1..consts::NT {
        // 1. Calculate velocity and position at current timestamp
        // ===
        // might also assume a linear electron density gradient?
        /*
        let dedendx = if meshx == mesh.nx - 1 { deden[[meshx, meshz]].0 } else {
            deden[[meshx, meshz]].0 + (
                (deden[[meshx+1, meshz]].0 - deden[[meshx, meshz]].0)
                / (mesh.get(meshx+1, meshz).x - mesh.get(meshx, meshz).x)
                * (x[tt-1] - mesh.get(meshx, meshz).x)
            )
        };

        // QUESTION: why is mydedendx computed but not mydedendz???
        // ANSWER: for the current test cases, because of the linear electron gradient,
        // mydedendz is constant across a row
        vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * dedendx * consts::DT;
        vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * deden[[meshx, meshz]].1 * consts::DT;*/

        // TODO interp
        vx -= consts::C_SPEED.powi(2) / (2.0 * consts::NCRIT) * deden[[meshx, meshy, meshz]].x * consts::DT;
        vy -= consts::C_SPEED.powi(2) / (2.0 * consts::NCRIT) * deden[[meshx, meshy, meshz]].y * consts::DT;
        vz -= consts::C_SPEED.powi(2) / (2.0 * consts::NCRIT) * deden[[meshx, meshy, meshz]].z * consts::DT;

        x.push(x[tt-1] + vx * consts::DT);
        y.push(y[tt-1] + vy * consts::DT);
        z.push(z[tt-1] + vz * consts::DT);

        // 2. update meshx, meshz (for vel calculation)
        /*(meshx, meshz) = mesh.get_mesh_coords_in_area(
            x[tt], z[tt],
            ( // min pt
                if meshx <= 1 { 0 } else { meshx - 1 },
                if meshz <= 1 { 0 } else { meshz - 1 },
            ),
            (min(mesh.nx - 1, meshx + 1), min(mesh.nz - 1, meshz + 1))
        );*/
        XYZ { x: meshx, y: meshy, z: meshz } = mesh.get_mesh_coords_in_area(
            XYZ { x: x[tt], y: y[tt], z: z[tt] },
            XYZ { x: max(meshx-1, 0), y: max(meshy-1, 0), z: max(meshz-1, 0) },
            XYZ { x: min(meshx+1, mesh.n.x-1), y: min(meshy+1, mesh.n.y-1), z: min(meshz+1, mesh.n.z-1) },
        );

        // 3. Stop the ray if it escapes the mesh
        if 
            x[tt] < mesh.min.x || x[tt] > mesh.max.x ||
            y[tt] < mesh.min.y || y[tt] > mesh.max.y ||
            z[tt] < mesh.min.z || z[tt] > mesh.max.z {
            break;
        }
    }

    (x, y, z)
}
