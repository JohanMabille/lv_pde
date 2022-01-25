#pragma once
#include <iostream>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include "Financial_PDE.hpp"

class BoundaryConditions
{
public:
	virtual ~BoundaryConditions() = default;

	Financial_PDE* pde_;

	const double& get_theta() const;
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const = 0;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const = 0;

protected:

	
	BoundaryConditions(Financial_PDE* pde, const double& theta);


private:

	double theta_;

};

class Boundaryx0: public BoundaryConditions 
{

public:

	Boundaryx0(Financial_PDE* pde, const double& theta);
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const override;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const override;

protected:

	
	

private:

	double gamma_x0(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double vega_x0(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double mu_x0(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
};

class BoundaryxN : public BoundaryConditions
{

public:

	BoundaryxN(Financial_PDE* pde, const double& theta);
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const override;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const override;

protected:

	

private:

	double gamma_xN(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double vega_xN(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double mu_xN(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
};

