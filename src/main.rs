//use cbet_rs::consts;
use cbet_rs::beam::*;
use cbet_rs::ray_trace;
use cbet_rs::cbet;

use std::time::SystemTime;

fn main() {
    println!("Creating and initializing mesh");
    println!("N rays: {}", cbet_rs::consts::NRAYS);
    //let mut m = cbet_rs::new_default_mesh();
    //m.init_points(consts::NCRIT, 2.0 * consts::NU_EI_C * 1e12/(consts::C_SPEED*1e4));
    //m.init_eden_machnum_3beam(consts::NCRIT);
    let m = cbet_rs::new_mesh();

    println!("Creating and initializing beams");
    let beam1 = Beam::beam1();
    //let beam3 = Beam::beam3();
    //let beam4 = Beam::beam4();
    //let mut beams = vec![beam1, beam3, beam4];
    let beam2 = Beam::beam2();
    let mut beams = vec![beam1, beam2];

    println!("Tracing rays");
    let now = SystemTime::now();
    ray_trace::ray_trace(&m, &mut beams);
    println!("\ttook {} seconds", now.elapsed().unwrap().as_secs_f64());

    println!("Doing CBET calculation");
    let now = SystemTime::now();
    cbet::init_crossings(&mut beams);
    cbet::cbet(&m, &mut beams);
    println!("\ttook {} seconds", now.elapsed().unwrap().as_secs_f64());

    println!("Saving to \"out.hdf5\"");
    cbet_rs::save_hdf5(&m, &beams, "out.hdf5").unwrap();
}
