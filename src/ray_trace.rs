use crate::beam::*;
use crate::mesh::*;
use crate::utils;
use crate::consts;
use std::cmp::min;
use std::sync::Mutex;
use rayon::prelude::*;

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
    let mut deden: Vec<(f64, f64)> = mesh.points.iter().enumerate().map(|(i, point)| {
        let x = i / mesh.nz;
        let z = i % mesh.nz;
        if x == mesh.nx - 1 || z == mesh.nz - 1 {
            (0.0, 0.0)
        } else {
            let dedendx = if x == 0 { 0.0 } else {
                (point.eden - mesh.get(x-1, z).eden)/(point.x-mesh.get(x-1, z).x)
            };
            let dedendz = if z == 0 { 0.0 } else {
                (point.eden - mesh.get(x, z-1).eden)/(point.z-mesh.get(x, z-1).z)
            };
            (dedendx, dedendz)
        }
    }).collect();

    // fill in outer points
    // ROW - i is the column in the row that's being modified
    for i in 0..mesh.nz {
        // first row: x = 0, fill in xes only
        // mesh.nz+i is point at (1, i)
        deden[i].0 = deden[mesh.nz+i].0;
        // fill in x = mesh.nx - 1
        // tuple implements copy so no need to worry
        deden[(mesh.nz*(mesh.nx-1))+i] = deden[(mesh.nz*(mesh.nx-2))+i];
    }
    // COL - i is the row in the column that's being modified
    for i in 0..mesh.nx {
        // first col: z = 0, fill in zs only
        deden[i*mesh.nz].1 = deden[i*mesh.nz+1].1;
        // fill in z = mesh.nz - 1
        deden[i*mesh.nz+(mesh.nz-1)] = deden[i*mesh.nz+(mesh.nz-2)];
    }

    // NOTE: wpe is a fn Mesh::wpe now instead of being precalculated
    // bad move? let me know.

    beams.iter_mut().for_each(|beam| {
        // Initialize the marked array
        beam.marked = (0..(mesh.nx*mesh.nz)).map(|_| Mutex::new(Vec::new())).collect();
        beam.rays.par_iter_mut().enumerate().for_each(|(i, ray)| {
            let child = launch_child_ray(ray, mesh, &deden);
            launch_parent_ray(ray, mesh, &deden, &child, &beam.marked, i);
        });
    });
}

/// Identical code for parent/child so I moved it out. Gets the k constant that is used to
/// calculate the initial velocities
fn get_k(mesh: &Mesh, x0: f64) -> f64 {
    // NOTE: it looks like this assumes a linear electron density gradient!
    // The Mesh::wpe function (currently commented out) finds the wpe without doing that,
    // but only for the point, not interpolated. Consider using that to find wpe instead?
    // Also, I changed pow(..., 2) to just multiplying it by itself
    let wpe_interp = f64::sqrt(
        utils::interp_closure(
            |i| { mesh.points[mesh.nz*i].eden }, mesh.nx, // eden.col(0) as y, ysize
            |i| { mesh.points[mesh.nz*i].x }, mesh.nx, // x.col(0) as x, xsize
            x0
        ) * 1e6 * consts::EC*consts::EC / (consts::ME*consts::E0)
    );
    // also, it looks like wpe_interp is only used here, and is sqrted only to be squared
    // again right after...
    // this is k
    f64::sqrt((consts::OMEGA*consts::OMEGA - wpe_interp*wpe_interp) / (consts::C_SPEED*consts::C_SPEED))
}

/// Launch the parent ray. Updates the marked vector, as well as the ray's crossings vector,
/// which it overwrites. More detailed comments for some of the code may be in
/// `launch_child_ray`, which I wrote first.
///
/// i is the ray index, because rays are stored in the marked array as indices (rather than
/// pointers, which might be what we could switch it to?)
fn launch_parent_ray(
    ray: &mut Ray, mesh: &Mesh, deden: &Vec<(f64, f64)>,
    (childx, childz): &(Vec<f64>, Vec<f64>), marked: &Vec<Mutex<Vec<(usize, usize)>>>, raynum: usize
) {
    // Maybe it is better to have all of these in one vector that stores a struct, i.e.
    // struct Timestamp { x, z, vx, vz, mx, mz }??
    // but I'll keep it like this for now, why not

    // renamed from myx, myz
    //let mut pos: Vec<(f64, f64)> = Vec::new();
    // renamed from myvx, myvz
    //let mut vel: Vec<(f64, f64)> = Vec::new();
    // looks like we only need these: + vx, vz later on
    let mut x = ray.x0;
    let mut z = ray.z0;

    // convert child coordinates into a set of distances (cumulative)
    // TODO: possibly move this code to launch_child_ray and do it in place?
    let mut distance = Vec::with_capacity(childx.len());
    distance.push(0.0);
    for i in 1..childx.len() {
        distance.push(distance[i-1] + f64::sqrt((childx[i] - childx[i-1]).powi(2) + (childz[i] - childz[i-1]).powi(2)));
    }

    let init_area = {
        let init_diff_x = ray.x0 - childx[0];
        let init_diff_z = ray.z0 - childz[0];
        let init_diff_mag = f64::sqrt(init_diff_x*init_diff_x + init_diff_z*init_diff_z);
        let init_proj_coeff = f64::abs(ray.kx0 * (init_diff_z/init_diff_mag) - ray.kz0 * (init_diff_x/init_diff_mag));
        init_diff_mag*init_proj_coeff
    };

    // renamed from thisx_0, thisz_0, as well as meshx, meshz
    let (mut meshx, mut meshz) = mesh.get_mesh_coords(ray.x0, ray.z0);

    let k = get_k(mesh, ray.x0);
    let knorm = f64::sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);

    let mut vx = consts::C_SPEED*consts::C_SPEED * ((ray.kx0 / knorm) * k) / consts::OMEGA;
    let mut vz = consts::C_SPEED*consts::C_SPEED * ((ray.kz0 / knorm) * k) / consts::OMEGA;

    let mut curr_dist = 0.0;
    for _ in 1..consts::NT {
        // 1. Calculate velocity and position at current timestamp
        // ===
        // **copied from launch_child_ray**
        let mydedendx = utils::interp_closure(
            |i| { deden[mesh.nz*i + (mesh.nz+1)/2].0 }, mesh.nx, //dedendx.col((nz+1)/2)
            |i| { mesh.points[mesh.nz*i].x }, mesh.nx, // x.col(0)
            x
        );

        // update pos, vel, equivalent to setting myx/z[tt] and myvx/z[tt]
        let prev_vx = vx;
        let prev_vz = vz;
        vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * mydedendx * consts::DT;
        vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * deden[meshx*mesh.nz + meshz].1 * consts::DT;
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
            // let crossx = utils::interp(&vec![prev_z, z], &vec![prev_x, x], currx);
            let crossx = prev_z + ((z - prev_z) / (x - prev_x) * (currx - prev_x)); // ***changed***
            let frac = (currx - prev_x) / (x - prev_x);
            assert!((frac >= 0.0) && (frac <= 1.0));

            if f64::abs(crossx - lastz) > 1.0e-12 {
                ray.crossings.push(Crossing {
                    x: currx,
                    z: crossx,
                    boxesx: meshx,
                    boxesz: meshz,
                    area_ratio: {
                        let extra = f64::sqrt((currx - prev_x).powi(2) + (crossx - prev_z).powi(2));
                        let childxp = utils::interp(childx, &distance, curr_dist + extra);
                        let childzp = utils::interp(childz, &distance, curr_dist + extra);

                        let diff_x = currx - childxp;
                        let diff_z = crossx - childzp;
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
                    w_mult: 1.0,
                });
                marked[meshx*mesh.nx + meshz].lock().unwrap().push((raynum, ray.crossings.len()-1));

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
            // let crossz = utils::interp(&vec![prev_x, x], &vec![prev_z, z], currz);
            let crossz = prev_x + ((x - prev_x) / (z - prev_z)) * (currz - prev_z);
            let frac = (currz - prev_z) / (z - prev_z);
            assert!((frac >= 0.0) && (frac <= 1.0));

            if f64::abs(crossz - lastx) > 1.0e-12 {
                ray.crossings.push(Crossing {
                    x: crossz,
                    z: currz,
                    boxesx: meshx,
                    boxesz: meshz,
                    area_ratio: {
                        let extra = f64::sqrt((crossz - prev_x).powi(2) + (currz - prev_z).powi(2));
                        let childxp = utils::interp(childx, &distance, curr_dist + extra);
                        let childzp = utils::interp(childz, &distance, curr_dist + extra);

                        let diff_x = crossz - childxp;
                        let diff_z = currz - childzp;
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
                    w_mult: 1.0,
                });
                marked[meshx*mesh.nx + meshz].lock().unwrap().push((raynum, ray.crossings.len()-1));
                is_cross_z = true;
                lastz = currz;
            }
        }
        // swap if out of order. then, calculate dkx, dkz, dkmag for prev. crossing(s)
        if is_cross_x && is_cross_z {
            let last_crossing_ind = ray.crossings.len()-1;
            if (x - prev_x) * (ray.crossings[last_crossing_ind].x - ray.crossings[last_crossing_ind-1].x) < 0.0 {
                ray.crossings.swap(last_crossing_ind, last_crossing_ind-1);
            }
            calc_dk(ray, last_crossing_ind-2);
            calc_dk(ray, last_crossing_ind-1);
        } else if (is_cross_x || is_cross_z) && ray.crossings.len() > 1 {
            calc_dk(ray, ray.crossings.len()-2);
        }
        //uray.push(uray[tt-1]);
        curr_dist += f64::sqrt((x - prev_x).powi(2) + (z - prev_z).powi(2));
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
    ray: &mut Ray, mesh: &Mesh, deden: &Vec<(f64, f64)>
) -> (Vec<f64>, Vec<f64>) {
    // renamed from myx, myz
    //let mut pos: Vec<(f64, f64)> = Vec::new();
    let mut x: Vec<f64> = Vec::new();
    let mut z: Vec<f64> = Vec::new();
    // renamed from myvx, myvz
    //let mut vel: Vec<(f64, f64)> = Vec::new();
    // Turns out, we only need to store current and previous velocity!
    x.push(ray.cx0);
    z.push(ray.cz0);

    // renamed from thisx_0, thisz_0
    let (mut meshx, mut meshz) = mesh.get_mesh_coords(ray.cx0, ray.cz0);

    let k = get_k(mesh, ray.cx0);
    let knorm = f64::sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);
    let mut vx = consts::C_SPEED*consts::C_SPEED * ((ray.kx0 / knorm) * k) / consts::OMEGA;
    let mut vz = consts::C_SPEED*consts::C_SPEED * ((ray.kz0 / knorm) * k) / consts::OMEGA;

    for tt in 1..consts::NT {
        // 1. Calculate velocity and position at current timestamp
        // ===
        // might also assume a linear electron density gradient?
        // TODO: move computing dedendx.col, x.col outside of for loop
        let mydedendx = utils::interp_closure(
            |i| { deden[mesh.nz*i + (mesh.nz+1)/2].0 }, mesh.nx, // dedendx.col((nz+1)/2)
            |i| { mesh.points[mesh.nz*i].x }, mesh.nx, // x.col(0)
            x[tt-1]
        );

        // QUESTION: why is mydedendx computed but not mydedendz???
        // ANSWER: for the current test cases, because of the linear electron gradient,
        // mydedendz is constant
        vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * mydedendx * consts::DT;
        vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0 * consts::NCRIT) * deden[meshx*mesh.nz + meshz].1 * consts::DT;

        x.push(x[tt-1] + vx * consts::DT);
        z.push(z[tt-1] + vz * consts::DT);

        // 2. update meshx, meshz (for vel calculation)
        (meshx, meshz) = mesh.get_mesh_coords_in_area(
            x[tt], z[tt],
            ( // min pt
                if meshx <= 1 { 0 } else { meshx - 1 },
                if meshz <= 1 { 0 } else { meshz - 1 },
            ),
            (min(mesh.nx - 1, meshx + 1), min(mesh.nz - 1, meshz + 1))
        );

        // 3. Stop the ray if it escapes the mesh
        if x[tt] < mesh.xmin || x[tt] > mesh.xmax || z[tt] < mesh.zmin || z[tt] > mesh.zmax {
            // Here, it outputs the trajectory to a file, TODO test???
            break;
        }
    }

    (x, z)
}
