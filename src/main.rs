use cbet_rs::consts;
use cbet_rs::beam::*;
use cbet_rs::ray_trace;
use cbet_rs::cbet;

fn main() {
    println!("Default mesh:");
    let mut m = cbet_rs::new_default_mesh();
    println!("Point at (50, 20): {:?}", m.get(50, 20));
    println!("Point at (55, 20): {:?}", m.get(55, 20));
    println!("Point at (50, 25): {:?}", m.get(50, 25));

    println!("Init eden, machnum:");
    m.init_eden_machnum(consts::NCRIT);
    println!("Point at (50, 20): {:?}", m.get(50, 20));
    println!("Point at (55, 20): {:?}", m.get(55, 20));
    println!("Point at (50, 25): {:?}", m.get(50, 25));

    let beam1 = Beam::beam1();
    let beam2 = Beam::beam2();
    let mut beams = vec![beam1, beam2];
    ray_trace::ray_trace(&m, &mut beams);

    //println!("{:?}", beams[0].rays[0]);

    cbet::cbet(&m, &mut beams);
}
