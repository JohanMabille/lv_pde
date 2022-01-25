#include "eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "pde_matrix_builder.hpp"


PDEMatrixBuilder::~PDEMatrixBuilder()
{
	delete boundaries_xN_;
	boundaries_xN_ = nullptr;
	delete boundaries_x0_;
	boundaries_x0_ = nullptr;
}

PDEMatrixBuilder::PDEMatrixBuilder(BoundaryConditions* boundaries_x0, BoundaryConditions* boundaries_xN) :
	boundaries_x0_(boundaries_x0), boundaries_xN_(boundaries_xN)
{

}

double PDEMatrixBuilder::alpha(const double& dx, const double& dt, const double& x, const double& t,
	const double& initial_vol, const double& initial_rate) const
{
	double alph = boundaries_x0_->pde_->coeff_c(x, t, initial_vol, initial_rate) - 
		2 * (boundaries_x0_->pde_->coeff_a(x, t, initial_vol, initial_rate)) / pow(dx, 2);
	return alph;
}

double PDEMatrixBuilder::beta(const double& dx, const double& dt, const double& x, const double& t,
	const double& initial_vol, const double& initial_rate) const
{
	double bet = (boundaries_x0_->pde_->coeff_a(x, t, initial_vol, initial_rate)) / pow(dx, 2) + 
		(boundaries_x0_->pde_->coeff_b(x, t, initial_vol, initial_rate)) / (2 * (dx));
	return bet;
}

double PDEMatrixBuilder::xi(const double& dx, const double& dt, const double& x, const double& t,
	const double& initial_vol, const double& initial_rate) const
{
	double result = (boundaries_x0_->pde_->coeff_a(x, t, initial_vol, initial_rate)) / pow(dx, 2) - 
		(boundaries_x0_->pde_->coeff_b(x, t, initial_vol, initial_rate)) / (2 * (dx));
	return result;
}

Eigen::MatrixXd PDEMatrixBuilder::comp_row_k(const double& dx, const double& dt, const double& x, const double& t,
	const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1, 3);
	res(0) = -dt * (1 - boundaries_x0_->get_theta()) * PDEMatrixBuilder::xi(dx, dt, x, t, initial_vol, initial_rate);
	res(1) = 1 - dt * (1 - boundaries_x0_->get_theta()) * PDEMatrixBuilder::alpha(dx, dt, x, t, initial_vol, initial_rate);
	res(2) = -dt * (1 - boundaries_x0_->get_theta()) * PDEMatrixBuilder::beta(dx, dt, x, t, initial_vol, initial_rate);
	return res;
}

Eigen::MatrixXd PDEMatrixBuilder::comp_row_ktilde(const double& dx, const double& dt, const double& x, const double& t,
	const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1, 3);
	res(0) = dt * boundaries_x0_->get_theta() * PDEMatrixBuilder::xi(dx, dt, x, t, initial_vol, initial_rate);
	res(1) = 1 + dt * boundaries_x0_->get_theta() * PDEMatrixBuilder::alpha(dx, dt, x, t, initial_vol, initial_rate);
	res(2) = dt * boundaries_x0_->get_theta() * PDEMatrixBuilder::beta(dx, dt, x, t, initial_vol, initial_rate);
	return res;
}

Eigen::MatrixXd PDEMatrixBuilder::comp_system_constant(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
	const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd vect = Eigen::MatrixXd::Zero(nb_spot_steps + 1, 1);
	double x;
	for (int i = 0; i <= nb_spot_steps; i++)
	{
		x = log_spot_prices(0, i);
		vect(i, 0) = dt * boundaries_x0_->pde_->coeff_d(x, t, initial_vol, initial_rate);
	} 
	return vect;
}

Eigen::MatrixXd PDEMatrixBuilder::comp_k(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
										const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(nb_spot_steps+1, nb_spot_steps+1); // null square matrix of dimension T+1
	//Filling the first row (boundary condition in x0):
	double x0 = log_spot_prices(0, 0);
	k.block<1, 3>(0, 0) = boundaries_x0_->coeff_fn1(dx, dt, x0, t, initial_vol, initial_rate);
	//filling the last row (boundary condition in xN):
	double xN = log_spot_prices(0, nb_spot_steps);
	k.block<1, 3>(nb_spot_steps, nb_spot_steps - 2) = boundaries_xN_->coeff_fn1(dx, dt, xN, t, initial_vol, initial_rate);

	double x;
	//loop to fill the intermediary rows:
	for (int i = 1; i < nb_spot_steps; i++)
	{
		x = log_spot_prices(0, i);
		k.block<1, 3>(i, i - 1) = PDEMatrixBuilder::comp_row_k(dx, dt, x, t, initial_vol, initial_rate);
	}
	return k;
}

Eigen::MatrixXd PDEMatrixBuilder::comp_ktilde(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
											const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd ktilde = Eigen::MatrixXd::Zero(nb_spot_steps +1, nb_spot_steps +1); // null square matrix of dimension T+1
	//Filling the first row (boundary condition in x0):
	double x0 = log_spot_prices(0, 0);
	ktilde.block<1, 3>(0, 0) = boundaries_x0_->coeff_fn(dx, dt, x0, t, initial_vol, initial_rate);
	//filling the last row (boundary condition in xN):
	double xN = log_spot_prices(0, nb_spot_steps);
	ktilde.block<1, 3>(nb_spot_steps, nb_spot_steps - 2) = boundaries_xN_->coeff_fn(dx, dt, xN, t, initial_vol, initial_rate);

	double x;
	//loop to fill the intermediary rows:
	for (int i = 1; i < nb_spot_steps; i++)
	{	
		x = log_spot_prices(0, i);
		ktilde.block<1, 3>(i, i - 1) = PDEMatrixBuilder::comp_row_ktilde(dx, dt, x, t, initial_vol, initial_rate);
	}
	return ktilde;
}