#pragma once

#include <iostream>
#include <cmath>

class VolatilityDiffusion
{
public:

	virtual ~VolatilityDiffusion() = default;
	// Pointeur vers une Option pour avoir la volatilité de départ
	//Option* Opt_;
	virtual double volatility_formula(const double& x, const double& t, const double& initial_vol) const = 0;

protected:

	VolatilityDiffusion() = default;
};

class ConstantVolatility : public VolatilityDiffusion
{

public:

	ConstantVolatility() = default;
	virtual double volatility_formula(const double& x, const double& t, const double& initial_vol) const override;

};