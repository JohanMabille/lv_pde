#pragma once

#include <iostream>
#include <cmath>

class RateDiffusion
{
public:

	virtual ~RateDiffusion() = default;
	// Pointeur vers une Option pour avoir la volatilité de départ
	//Option* Opt_;
	virtual double rate_formula(const double& x, const double& t, const double& initial_rate) const = 0;

protected:

	RateDiffusion() = default;
};

class ConstantRate : public RateDiffusion
{

public:

	ConstantRate() = default;
	virtual double rate_formula(const double& x, const double& t, const double& initial_rate) const override;

};