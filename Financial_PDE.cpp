#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <math.h> 
#include "Financial_PDE.hpp"



Financial_PDE::Financial_PDE(VolatilityDiffusion* vol, RateDiffusion* rate) : vol_(vol), rate_(rate)
{

}

BS_PDE::BS_PDE(VolatilityDiffusion* vol, RateDiffusion* rate) : Financial_PDE(vol, rate)
{
}

double BS_PDE::coeff_a(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const {
	double sigma = vol_->volatility_formula(x, t, initial_vol);
	return -0.5 * pow(sigma, 2);
}
double BS_PDE::coeff_b(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const {
	double sigma = vol_->volatility_formula(x, t, initial_vol);
	double r = rate_->rate_formula(x, t, initial_rate);
	return 0.5 * pow(sigma,2) - r;
}
double BS_PDE::coeff_c(const double& x, const double& t, const double& initial_vol, const double& initial_rate)const {
	double r = rate_->rate_formula(x, t, initial_rate);
	return r;
}
double BS_PDE::coeff_d(const double& x, const double& t, const double& initial_vol, const double& initial_rate)const {
	return 0;
}