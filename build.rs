fn main() {
    cxx_build::bridge("src/cbet.rs")
        .cuda(true)
        .file("src/cpp_cbet.cu")
        .compile("cbet_rs");

    println!("cargo:rerun-if-changed=src/cbet.rs");
    println!("cargo:rerun-if-changed=src/cpp_cbet.cc");
    println!("cargo:rerun-if-changed=include/cpp_cbet.cc");
}
