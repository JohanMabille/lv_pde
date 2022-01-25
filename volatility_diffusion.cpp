#include "volatility_diffusion.hpp"
#include <iostream>


double ConstantVolatility::volatility_formula(const double& x, const double& t, const double& initial_vol) const
{
	double vol = initial_vol;
	return vol;
}
