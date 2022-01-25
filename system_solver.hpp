#pragma once

#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"
#include "pde_matrix_builder.hpp"

class SystemSolver 
{
public:

	~SystemSolver();
	SystemSolver(PDEMatrixBuilder* mat_build);
	PDEMatrixBuilder* mat_build_;
	Eigen::MatrixXd get_fn(Eigen::MatrixXd& fn1, const int& nb_spot_steps, const double& dx, const double& dt, const double& t, const double& t1,
							const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const;
	Eigen::MatrixXd tridiagonal_inversion(const Eigen::MatrixXd& A) const;

private:

	Eigen::MatrixXd inverse_ktilde(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
									const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const;
	Eigen::MatrixXd thomas_algo(const Eigen::MatrixXd& A, const Eigen::MatrixXd& d) const;
	

};