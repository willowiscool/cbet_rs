#include "cbet_rs/include/cpp_cbet.h"
#include "cbet_rs/src/cbet.rs.h"
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

// TODO: MOVE TO USE RUST VALUES!!!
// [brc]n = beam/ray/crossing number
// n[brc] = number of beams/rays/crossings

const double CONVERGE = 1e-7;
const double MAX_INCR = 0.2;
const double CBETCONVERGENCE = 0.9990;

void cpp_cbet(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, size_t nb, size_t nr, size_t nc) {
	printf("Running cbet in cpp\n");
	// create w_mult_values struct:
	double* w_mult_values = new double[nb*nr*nc];

	double currmax = MAX_INCR;
	for (size_t i = 1; i <= 500; i++) {
		get_cbet_gain(cbet_crossings, cbet_crosses, w_mult_values, nb, nr, nc);

		double updateconv = update_intensities(cbet_crossings, w_mult_values, nb, nr, nc, 0.0, currmax);
		if (updateconv <= CONVERGE) {
			break;
		}

		double currmaxa = MAX_INCR*pow(CBETCONVERGENCE, i);
		double currmaxb = CBETCONVERGENCE*updateconv;
		currmax = std::min(currmaxa, currmaxb);
	}
}

void get_cbet_gain(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, double* w_mult_values, size_t nb, size_t nr, size_t nc) {
	for (size_t bn = 0; bn < nb; bn++) {
		for (size_t rn = 0; rn < nr; rn++) {
			for (size_t cn = 0; cn < nc; cn++) {
				CbetCrossing crossing = cbet_crossings[(((bn*nr)+rn)*nc)+cn];
				if (crossing.intensity == 0.0) break;
				double cbet_sum = 0.0;
				for (size_t crossn = 0; crossn < nb; crossn++) {
					// xs = crosses
					CbetCrosses xs = cbet_crosses[(((((bn*nr)+rn)*nc)+cn)*nb)+crossn];
					if (xs.coupling_mult == 0.0) break;
					double other_intensity1 = cbet_crossings[(((xs.b_num*nr)+xs.r_num)*nc)+xs.c_num].intensity;
					double other_intensity2 = cbet_crossings[(((xs.b_num*nr)+xs.r_num)*nc)+xs.c_num_next].intensity;
					double avg_intensity = (other_intensity1+other_intensity2)/2.0;
					cbet_sum += avg_intensity*xs.coupling_mult;
				}
				w_mult_values[(((bn*nr)+rn)*nc)+cn] = exp(-1.0*cbet_sum) * crossing.absorption_coeff;
			}
		}
	}
}
double update_intensities(CbetCrossing* cbet_crossings, double* w_mult_values, size_t nb, size_t nr, size_t nc, double conv_max, double curr_max) {
	double curr_conv_max = conv_max;
	for (size_t bn = 0; bn < nb; bn++) {
		for (size_t rn = 0; rn < nr; rn++) {
			double i0 = cbet_crossings[(((bn*nr)+rn)*nc)+0].intensity;
			double mult_acc = 1.0;
			for (size_t cn = 0; cn < nc; cn++) {
				// limit energy
				CbetCrossing crossing = cbet_crossings[(((bn*nr)+rn)*nc)+cn];
				if (crossing.intensity == 0.0) break;
				double i_curr = i0*mult_acc;
				double fractional_change = abs(i_curr-crossing.intensity)/crossing.intensity;
				curr_conv_max = std::max(fractional_change, curr_conv_max);
				if (fractional_change > curr_max) {
					int sign = i_curr - crossing.intensity > 0.0 ? 1.0 : -1.0;
					double correction = 1.0 + curr_max*sign;
					i_curr = crossing.intensity*correction;
				}
				// end limit energy
				mult_acc *= w_mult_values[(((bn*nr)+rn)*nc)+cn];
				crossing.intensity = i_curr;
				cbet_crossings[(((bn*nr)+rn)*nc)+cn] = crossing;
			}
		}
	}
	return curr_conv_max;
}
