#include "cbet_rs/include/cpp_cbet.cuh"
#include "cbet_rs/src/cbet.rs.h"
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

#define CEIL_DIV(a, b) ((a+b-1)/b)

// TODO: MOVE TO USE RUST VALUES!!!
// [brc]n = beam/ray/crossing number
// n[brc] = number of beams/rays/crossings

const double CONVERGE = 1e-7;
const double MAX_INCR = 0.2;
const double CBETCONVERGENCE = 0.9990;

// headers that can't go in .h file bc compiler reasons
__global__ void get_cbet_gain(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, double* w_mult_values, size_t nb, size_t nr, size_t nc);
__global__ void update_intensities(CbetCrossing* cbet_crossings, double* w_mult_values, size_t nb, size_t nr, size_t nc, double* conv_max, double curr_max);

void cpp_cbet(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, size_t nb, size_t nr, size_t nc) {
	printf("Running cbet in cpp\n");
	// create w_mult_values struct:
	//double* w_mult_values = new double[nb*nr*nc];

	CbetCrossing* cuda_cbet_crossings;
	CbetCrosses* cuda_cbet_crosses;
	double* cuda_w_mult_values;
	double* cuda_conv_max;
	double conv_max = 0.0;
	cudaMalloc(&cuda_cbet_crossings, nb*nr*nc * sizeof(CbetCrossing));
	cudaMalloc(&cuda_cbet_crosses, nb*nr*nc*nb * sizeof(CbetCrosses));
	cudaMalloc(&cuda_w_mult_values, nb*nr*nc * sizeof(double));
	cudaMalloc(&cuda_conv_max, sizeof(double));
	cudaMemcpy(cuda_cbet_crossings, cbet_crossings, nb*nr*nc * sizeof(CbetCrossing), cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_cbet_crosses, cbet_crosses, nb*nr*nc*nb * sizeof(CbetCrosses), cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_conv_max, &conv_max, sizeof(double), cudaMemcpyHostToDevice);

	// let's say each thread processes a single ray
	// there are nb*nr rays
	dim3 threads_per_block(32, 32);
	int blocks = CEIL_DIV(nb*nr, 1024);

	double currmax = MAX_INCR;
	for (size_t i = 1; i <= 500; i++) {
		conv_max = 0.0;
		cudaMemcpy(cuda_conv_max, &conv_max, sizeof(double), cudaMemcpyHostToDevice);

		get_cbet_gain<<<blocks, threads_per_block>>>(cuda_cbet_crossings, cuda_cbet_crosses, cuda_w_mult_values, nb, nr, nc);
		update_intensities<<<blocks, threads_per_block>>>(cuda_cbet_crossings, cuda_w_mult_values, nb, nr, nc, cuda_conv_max, currmax);

		cudaMemcpy(&conv_max, cuda_conv_max, sizeof(double), cudaMemcpyDeviceToHost);
		if (conv_max <= CONVERGE) {
			break;
		}

		double currmaxa = MAX_INCR*pow(CBETCONVERGENCE, i);
		double currmaxb = CBETCONVERGENCE*conv_max;
		currmax = std::min(currmaxa, currmaxb);
	}

	cudaMemcpy(cbet_crossings, cuda_cbet_crossings, nb*nr*nc * sizeof(CbetCrossing), cudaMemcpyDeviceToHost);
}

__global__ void get_cbet_gain(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, double* w_mult_values, size_t nb, size_t nr, size_t nc) {
	// ray index = beam number * nr + ray number
	size_t ray_index = (((blockIdx.x*1024)+threadIdx.y)*32)+threadIdx.x;
	if (ray_index > nb*nr) return;
	for (size_t cn = 0; cn < nc; cn++) {
		// should take address rather than value?
		CbetCrossing crossing = cbet_crossings[(ray_index*nc)+cn];
		if (crossing.intensity == 0.0) break;
		double cbet_sum = 0.0;
		for (size_t crossn = 0; crossn < nb; crossn++) {
			// xs = crosses
			CbetCrosses xs = cbet_crosses[(((ray_index*nc)+cn)*nb)+crossn];
			if (xs.coupling_mult == 0.0) break;
			double other_intensity1 = cbet_crossings[(((xs.b_num*nr)+xs.r_num)*nc)+xs.c_num].intensity;
			double other_intensity2 = cbet_crossings[(((xs.b_num*nr)+xs.r_num)*nc)+xs.c_num_next].intensity;
			double avg_intensity = (other_intensity1+other_intensity2)/2.0;
			cbet_sum += avg_intensity*xs.coupling_mult;
		}
		w_mult_values[(ray_index*nc)+cn] = exp(-1.0*cbet_sum) * crossing.absorption_coeff;
	}
}

// changes conv max to a pointer!!! writes into it!!!
__global__ void update_intensities(CbetCrossing* cbet_crossings, double* w_mult_values, size_t nb, size_t nr, size_t nc, double* conv_max, double curr_max) {
	size_t ray_index = (((blockIdx.x*1024)+threadIdx.y)*32)+threadIdx.x;
	if (ray_index > nb*nr) return;
	double i0 = cbet_crossings[(ray_index*nc)+0].intensity;
	double mult_acc = 1.0;
	for (size_t cn = 0; cn < nc; cn++) {
		CbetCrossing crossing = cbet_crossings[(ray_index*nc)+cn];
		if (crossing.intensity == 0.0) break;
		double i_curr = i0*mult_acc;
		double fractional_change = abs(i_curr-crossing.intensity)/crossing.intensity;
		*conv_max = max(fractional_change, *conv_max);
		if (fractional_change > curr_max) {
			int sign = i_curr - crossing.intensity > 0.0 ? 1.0 : -1.0;
			double correction = 1.0 + curr_max*sign;
			i_curr = crossing.intensity*correction;
		}
		mult_acc *= w_mult_values[(ray_index*nc)+cn];
		crossing.intensity = i_curr;
		cbet_crossings[(ray_index*nc)+cn] = crossing;
	}
}
