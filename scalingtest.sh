#!/usr/bin/bash

echo strong scaling test!
cargo build -r
threads=$(($(lscpu -e | wc -l)-1))
export SCALING_FACTOR=1
for i in $(seq 1 $threads)
do
	echo "$i thread(s)"
	time RAYON_NUM_THREADS="$i" target/release/cbet_rs > /dev/null
	time RAYON_NUM_THREADS="$i" target/release/cbet_rs > /dev/null
	time RAYON_NUM_THREADS="$i" target/release/cbet_rs > /dev/null
done
