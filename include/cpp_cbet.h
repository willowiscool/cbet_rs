#pragma once
#include <cstddef>

struct CbetCrosses;
struct CbetCrossing;

void cpp_cbet(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, size_t nb, size_t nr, size_t nc);

void get_cbet_gain(CbetCrossing* cbet_crossings, CbetCrosses* cbet_crosses, double* w_mult_values, size_t nb, size_t nr, size_t nc);
double update_intensities(CbetCrossing* cbet_crossings, double* w_mult_values, size_t nb, size_t nr, size_t nc, double conv_max, double curr_max);

