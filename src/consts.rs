// Program constants copied from def.h
// ===
// Mesh constants
// usize type is used for integers because they are often the lengths of vectors, which are usize
pub const NX: usize = 1001;
pub const XMIN: f64 = -8.0e-4;
pub const XMAX: f64 = 8.0e-4;

pub const NZ: usize = 1001;
pub const ZMIN: f64 = -8.0e-4;
pub const ZMAX: f64 = 8.0e-4;

// Beam constants
pub const RAYS_PER_ZONE: usize = 4;
pub const BEAM_MAX_Z: f64 = 3.0e-4;
pub const BEAM_MIN_Z: f64 = -3.0e-4;
// replaced DZ with the derivation for DZ
// the as usize operation is not necessarily safe but it's fine I guess
pub const NRAYS: usize = (RAYS_PER_ZONE as f64 * (BEAM_MAX_Z-BEAM_MIN_Z)/((ZMAX-ZMIN)/(NZ as f64-1.0))) as usize;
pub const OFFSET1: f64 = 0.1e-4;
pub const OFFSET2: f64 = 0.2e-4;
//pub const OFFSET2: f64 = 0.1e-4;
pub const DIR1: [f64; 2] = [1.0, 0.0];
pub const DIR2: [f64; 2] = [0.0, 1.0];
pub const CHILD_OFFSET: f64 = 0.1e-4;
// The following definitions are all for uray_mult
// comments copied from c++ impl.
pub const SIGMA: f64 = 2.0e-4;
pub const INTENSITY: f64 = 1e17; // intensity of the beam in W/cm^2
pub const COURANT_MULT: f64 = 0.1; // 0.37 // 0.25 // 0.36 // 0.22;
pub const URAY_MULT: f64 = INTENSITY*COURANT_MULT*(1.0/RAYS_PER_ZONE as f64);
// pub const URAY_MULT1: f64 = 5.0*INTENSITY*COURANT_MULT*(1.0/RAYS_PER_ZONE as f64);
pub const URAY_MULT2: f64 = 2.5*INTENSITY*COURANT_MULT*(1.0/RAYS_PER_ZONE as f64);

// Ray tracing constants
// https://stackoverflow.com/questions/53619695/calculating-maximum-value-of-a-set-of-constant-expressions-at-compile-time
const fn max_usize(a: usize, b: usize) -> usize {
    [a, b][(a < b) as usize]
}
/* Const functions do not support floating point arithmetic at compile time
const fn max_f64(a: f64, b: f64) -> f64 {
    [a, b][(a < b) as usize]
}*/
pub const NT: usize = ((1.0 / COURANT_MULT) * max_usize(NX, NZ) as f64 * 2.0) as usize;
// replaced DX, DZ with their derivations
//pub const DT: f64 = COURANT_MULT*max_f64((XMAX-XMIN)/(NX as f64-1.0), (ZMAX-ZMIN)/(NZ as f64-1.0))/C_SPEED;
// because I can't use max_f64, turns out we're defining dx and dz consts anyway, for now.
const DX: f64 = (XMAX-XMIN)/(NX as f64-1.0);
const DZ: f64 = (ZMAX-ZMIN)/(NZ as f64-1.0);
pub const DT: f64 = COURANT_MULT*([DX, DZ][(DX < DZ) as usize])/C_SPEED;

// CBET constants
pub const MAX_INCR: f64 = 0.2;
pub const CONVERGE: f64 = 1e-7;
pub const CBETCONVERGENCE: f64 = 0.9990;

// ===
// scientific consts copied from def.h, comments also copied haha
// ===
// TODO: consider adding a C_SPEED^2 const
pub const C_SPEED: f64 = 29979245800.0; // speed of light in cm/s
pub const E0: f64 = 8.85418782e-12; // permittivity of free space in m^-3 kg^-1 s^4 A^2
pub const ME: f64 = 9.10938356e-31; // electron mass in kg
pub const EC: f64 = 1.60217662e-19; // electron charge in C

pub const LAMBDA: f64 = 1.053e-4/3.0; // wavelength of light, in cm. This is frequency-tripled
                                      // "3w" or "blue" (UV) light
pub const FREQ: f64 = C_SPEED/LAMBDA; // frequency of light, in Hz
pub const OMEGA: f64 = 2.0*std::f64::consts::PI*FREQ; // frequency of light, in rad/s
// the critical density occurs when omega = omega_p,e
// NOTE: replaced pow(omega,2) with omega*omega and pow(ec, 2) with ec*ec
pub const NCRIT: f64 = 1e-6*(OMEGA*OMEGA*ME*E0/(EC*EC));

// More CBET constants
pub const ESTAT: f64 = 4.80320427e-10; // electron charge in statC
pub const ME_G: f64 = 9.10938356e-28; // electron mass in g
pub const KB: f64 = 1.3806485279e-16; // Boltzmann constant in erg/K
pub const NORM: f64 = 1e14;
pub const CBET_CONST: f64 = (8.0*std::f64::consts::PI*1e7*NORM/C_SPEED)*(ESTAT*ESTAT/(4.0*ME_G*C_SPEED*KB*1.1605e7))*1e-4;

pub const Z: f64 = 3.1; // ionization state
pub const TE_EV: f64 = 2.0e3; // Te_eV; the comment on Te is "Temperature of electron in K"
pub const TI_EV: f64 = 1.0e3; // Ti_eV; the comment on Ti is "Temperature of ion in K"
pub const MI_KG: f64 = 10230.0*ME; // Mass of ion in kg
// uh oh not a real const because of sqrt
pub static mut CS: f64 = 0.0;
pub const IAW: f64 = 0.542940629585429; // ion-acoustic wave energy-damping rate (nu_ia/omega_s)!!

pub fn init_consts() {
    unsafe {
        CS = 1e2*f64::sqrt(EC*(Z*TE_EV+3.0*TI_EV)/MI_KG);
    }
}
