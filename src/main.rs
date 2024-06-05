use cbet_rs::consts;
use cbet_rs::beam::*;
use cbet_rs::ray_trace;
use cbet_rs::cbet;

use std::time::SystemTime;

fn main() {
    println!("Creating and initializing mesh");
    let mut m = cbet_rs::new_default_mesh();
    m.init_eden_machnum(consts::NCRIT);

    println!("Creating and initializing beams");
    let beam1 = Beam::beam1();
    let beam2 = Beam::beam2();
    let mut beams = vec![beam1, beam2];

    println!("Tracing rays");
    let now = SystemTime::now();
    ray_trace::ray_trace(&m, &mut beams);
    println!("\ttook {} seconds", now.elapsed().unwrap().as_secs_f64());

    println!("Doing CBET calculation");
    let now = SystemTime::now();
    cbet::init_crossings(&mut beams, consts::INTENSITY);
    cbet::cbet(&m, &mut beams);
    println!("\ttook {} seconds", now.elapsed().unwrap().as_secs_f64());

    println!("Saving to \"out.hdf5\"");
    cbet_rs::save_hdf5(&m, &beams, "out.hdf5").unwrap();
}
