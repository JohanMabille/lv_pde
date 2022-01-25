#pragma once
#include "eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <algorithm>
#include "boundary_conditions.hpp"


class PDEMatrixBuilder
{

public:

	~PDEMatrixBuilder();
	PDEMatrixBuilder(BoundaryConditions* boundaries_x0, BoundaryConditions* boundaries_xN);
	Eigen::MatrixXd comp_k(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
							const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const;
	Eigen::MatrixXd comp_ktilde(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
								const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const;
	BoundaryConditions* boundaries_x0_;
	BoundaryConditions* boundaries_xN_;

	Eigen::MatrixXd comp_system_constant(const int& nb_spot_steps, const double& dx, const double& dt, 
							const double& t, const Eigen::MatrixXd& log_spot_prices,
							const double& initial_vol, const double& initial_rate) const; //compute the constant vector in the system

protected:

	
	
private:

	double alpha(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double beta(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	double xi(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;

	
	Eigen::MatrixXd comp_row_ktilde(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	Eigen::MatrixXd comp_row_k(const double& dx, const double& dt, const double& x, const double& t,
		const double& initial_vol, const double& initial_rate) const;
	

};

