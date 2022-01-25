#include "eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "pde_matrix_builder.hpp"
#include "system_solver.hpp"
#include <chrono>

SystemSolver::~SystemSolver()
{
	delete mat_build_;
	mat_build_ = nullptr;
}

SystemSolver::SystemSolver(PDEMatrixBuilder* mat_build) : mat_build_(mat_build)
{

}

Eigen::MatrixXd SystemSolver::inverse_ktilde(const int& nb_spot_steps, const double& dx, const double& dt, const double& t,
											const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const
{
	Eigen::MatrixXd ktilde = mat_build_->comp_ktilde(nb_spot_steps, dx, dt, t, log_spot_prices, initial_vol, initial_rate);
	Eigen::MatrixXd res = ktilde.inverse();
	return res;
}

Eigen::MatrixXd SystemSolver::get_fn(Eigen::MatrixXd& fn1, const int& nb_spot_steps, const double& dx, const double& dt, const double& t, const double& t1,
									const Eigen::MatrixXd& log_spot_prices, const double& initial_vol, const double& initial_rate) const
{	
	Eigen::MatrixXd k = mat_build_->comp_k(nb_spot_steps, dx, dt, t1, log_spot_prices, initial_vol, initial_rate);
	Eigen::MatrixXd ktilde = mat_build_->comp_ktilde(nb_spot_steps, dx, dt, t, log_spot_prices, initial_vol, initial_rate);
	Eigen::MatrixXd fn = Eigen::MatrixXd::Zero(nb_spot_steps + 1, 1);
	Eigen::MatrixXd ktilde_inv = SystemSolver::tridiagonal_inversion(ktilde);
	/*
	fn = SystemSolver::inverse_ktilde(nb_spot_steps, dx, dt, t, log_spot_prices, initial_vol, initial_rate) * k * fn1 -
						dt * k * mat_build_->comp_system_constant(nb_spot_steps, dx, dt, t, log_spot_prices, initial_vol, initial_rate);
	*/
	fn = ktilde_inv * k * fn1 -
		dt * ktilde_inv * (mat_build_->boundaries_x0_->get_theta() * mat_build_->comp_system_constant(nb_spot_steps, dx, dt, t, log_spot_prices, initial_vol, initial_rate) + 
		(1 - mat_build_->boundaries_x0_->get_theta()) * mat_build_->comp_system_constant(nb_spot_steps, dx, dt, t1, log_spot_prices, initial_vol, initial_rate));
	
	return fn;
}

Eigen::MatrixXd SystemSolver::thomas_algo(const Eigen::MatrixXd& A, const Eigen::MatrixXd& d) const
{
	int dimension = A.rows();

	Eigen::MatrixXd a = Eigen::MatrixXd::Zero(dimension, 1);
	Eigen::MatrixXd b = Eigen::MatrixXd::Zero(dimension, 1);
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(dimension, 1);
	Eigen::MatrixXd c_prime = Eigen::MatrixXd::Zero(dimension, 1);
	Eigen::MatrixXd d_prime = Eigen::MatrixXd::Zero(dimension, 1);

	//filling the vectors:
	a(0) = 0;
	a(dimension - 1) = A(dimension - 1, dimension - 2);
	b(0) = A(0, 0);
	b(dimension - 1) = A(dimension - 1, dimension - 1);
	c(0) = A(0, 1);
	c(dimension - 1) = 0;

	c_prime(0) = c(0) / b(0);
	d_prime(0) = d(0) / b(0);

	double denominator;

	for (int i = 1; i < dimension - 1; i++)
	{
		a(i) = A(i, i - 1);
		b(i) = A(i, i);
		c(i) = A(i, i + 1);
		denominator = (b(i) - a(i) * c_prime(i - 1));
		c_prime(i) = c(i) / denominator;
		d_prime(i) = (d(i) - a(i) * d_prime(i - 1)) / denominator;
	}
	denominator = (b(dimension-1) - a(dimension - 1) * c_prime(dimension - 2));
	c_prime(dimension - 1) = c(dimension - 1) / denominator;
	d_prime(dimension - 1) = (d(dimension - 1) - a(dimension - 1) * d_prime(dimension - 2)) / denominator;
	

	Eigen::MatrixXd sol = Eigen::MatrixXd::Zero(dimension, 1);
	sol(dimension - 1) = d_prime(dimension - 1);

	for (int i = 1; i < dimension - 1; i++)
	{
		sol(dimension - 1 - i) = d_prime(dimension - 1 - i) - c_prime(dimension - 1 - i) * sol(dimension - i);
	}

	return sol;
}

Eigen::MatrixXd SystemSolver::tridiagonal_inversion(const Eigen::MatrixXd& A) const
{
	int dimension = A.rows();
	Eigen::MatrixXd inverse_mat = Eigen::MatrixXd::Zero(dimension, dimension);
	for (int i = 0; i < dimension; i++)
	{
		Eigen::MatrixXd d = Eigen::MatrixXd::Zero(dimension, 1);
		d(i) = 1;
		inverse_mat.col(i) = SystemSolver::thomas_algo(A, d);
	}

	return inverse_mat;

}