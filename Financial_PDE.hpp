#pragma once
#ifndef FINANCIAL_PDE_HPP
#define FINANCIAL_PDE_HPP

#include "volatility_diffusion.hpp"
#include "rate_diffusion.hpp"


class Financial_PDE {

public:

	virtual ~Financial_PDE() = default;

	Financial_PDE(VolatilityDiffusion* vol, RateDiffusion* rate);
	// Or "Option* Opt;" ??

	//=== coefficients
	virtual double coeff_a(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const = 0; // ou (double t, double x)
	virtual double coeff_b(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const = 0;
	virtual double coeff_c(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const = 0;
	virtual double coeff_d(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const = 0;
	VolatilityDiffusion* vol_;
	RateDiffusion* rate_;

protected:
	

private:



};

class BS_PDE :public Financial_PDE {

public:

	BS_PDE(VolatilityDiffusion* vol, RateDiffusion* rate);


	//=== BS coefficients


	double coeff_a(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const; // ou (double t, double x)
	double coeff_b(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const;
	double coeff_c(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const;
	double coeff_d(const double& x, const double& t, const double& initial_vol, const double& initial_rate) const;

protected:

private:


};

#endif
