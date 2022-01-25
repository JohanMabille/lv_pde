#include "eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "system_solver.hpp"
#include "mesh.hpp"
// Solving execution time problems.
#include <chrono>

Mesh::Mesh(SystemSolver* syst_solv, Payoff* payoff, const int& nb_spot_steps, const int& nb_time_steps, const double& initial_spot,
	const double& K, const double& initial_vol, const double& T, const double& initial_rate, const bool& constant_coeffs, const double& eps) :
	syst_solv_(syst_solv), payoff_(payoff), nb_spot_steps_(nb_spot_steps), nb_time_steps_(nb_time_steps), initial_spot_(initial_spot),
	K_(K), initial_vol_(initial_vol), T_(T), initial_rate_(initial_rate), constant_coeffs_(constant_coeffs), eps_(eps)
{

}

double Mesh::get_K() const
{
	return K_;
}
double Mesh::get_vol() const
{
	return initial_vol_;
}
double Mesh::get_T() const
{
	return T_;
}
double Mesh::get_rate() const
{
	return initial_rate_;
}

const int& Mesh::get_nb_time_steps() const
{
	return nb_time_steps_;
}

const int& Mesh::get_nb_spot_steps() const
{
	return nb_spot_steps_;
}

const double& Mesh::get_initial_spot() const
{
	return initial_spot_;
}

double Mesh::comp_dt() const
{
	double maturity = T_;
	double dt = maturity / nb_time_steps_;
	return dt;
}

double Mesh::comp_dx() const
{
	double std = initial_vol_ * pow(T_, 0.5);
	double width = 10*std;
	double dx = width / nb_spot_steps_;
	return dx;
}

Eigen::MatrixXd Mesh::log_spot_prices() const
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1, nb_spot_steps_ + 1);
	double std = initial_vol_ * pow(T_, 0.5);
	double initial_log_spot = log(initial_spot_);
	double x0 = initial_log_spot - 5 * std;
	res(0, 0) = x0;
	double dx = Mesh::comp_dx();
	for (int i = 1; i <= nb_spot_steps_; i++)
	{
		res(0, i) = x0 + i * dx;
	}
	return res;
}


Eigen::MatrixXd Mesh::mesh_maturity() const
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1, nb_spot_steps_ + 1);
	Eigen::MatrixXd log_spots = Mesh::log_spot_prices();
	for (int i = 0; i <= nb_spot_steps_; i++)
	{
		double log_spot = log_spots(0, i);
		res(0, i) = payoff_->operator()(exp(log_spot));
	}
	return res;
}

Eigen::MatrixXd Mesh::comp_mesh(const Eigen::MatrixXd& log_spots, const bool bumped) const
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(nb_time_steps_+1, nb_spot_steps_ + 1);
	//res.block<1, nb_time_steps_ + 1>(0, 0) = Mesh::mesh_maturity();
	res.row(0) = Mesh::mesh_maturity();
	Eigen::MatrixXd fn1 = Eigen::MatrixXd::Zero(nb_spot_steps_ + 1, 1);
	Eigen::MatrixXd solved = Eigen::MatrixXd::Zero(nb_spot_steps_ + 1, 1);
	double dx = Mesh::comp_dx();
	double dt = Mesh::comp_dt();

	double initial_vol_used;
	double t;
	double t1;

	if (bumped == true)
	{
		initial_vol_used = initial_vol_ + eps_;
	}
	else
	{
		initial_vol_used = initial_vol_;
	}

	if (constant_coeffs_ == true)
	{	
		Eigen::MatrixXd ktilde = syst_solv_->mat_build_->comp_ktilde(nb_spot_steps_, dx, dt, t, log_spots, initial_vol_used, initial_rate_);
		//Eigen::MatrixXd ktilde_inv = syst_solv_->tridiagonal_inversion(ktilde);
		Eigen::MatrixXd ktilde_inv = ktilde.inverse();
		Eigen::MatrixXd k = syst_solv_->mat_build_->comp_k(nb_spot_steps_, dx, dt, t, log_spots, initial_vol_used, initial_rate_);
		Eigen::MatrixXd system_constant = dt * ktilde_inv * syst_solv_->mat_build_->comp_system_constant(nb_spot_steps_, dx, dt, t, log_spots, initial_vol_used, initial_rate_);

		Eigen::MatrixXd ktildeinv_k_product = ktilde_inv * k;

		for (int i = 1; i <= nb_time_steps_; i++)
		{
			
			fn1 = res.row(i - 1).transpose();
			
			solved = ktildeinv_k_product * fn1 - system_constant;
			res.row(i) = solved.transpose();
		}
	}
	else
	{
		for (int i = 1; i <= nb_time_steps_; i++)
		{
			t1 = T_ - (i-1) * dt;
			t = T_ - i * dt;
			fn1 = res.row(i - 1).transpose();
			solved = syst_solv_->get_fn(fn1, nb_spot_steps_, dx, dt, t, t1, log_spots, initial_vol_used, initial_rate_);
			res.row(i) = solved.transpose();
		}
	}
	return res;
}

double Mesh::linear_interp_around_spot(const Eigen::MatrixXd& log_spots, const double& y_inf, const double& y_sup) const
{
	int position_inf = (nb_spot_steps_ - 1) / 2;
	int position_sup = (nb_spot_steps_ - 1) / 2 + 1;
	double x_inf = log_spots(position_inf);
	double x_sup = log_spots(position_sup);
	double res = ((x_sup - log(initial_spot_)) * y_inf + (log(initial_spot_) - x_inf) * y_sup) / (x_sup - x_inf);
	return res;
}

double Mesh::extract_price(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const
{
	double price;
	if (nb_spot_steps_ % 2 == 0)
	{
		price = full_mesh(nb_time_steps_, nb_spot_steps_ / 2);
	}
	else
	{
		int position_inf = (nb_spot_steps_ - 1) / 2;
		int position_sup = (nb_spot_steps_ - 1) / 2 + 1;
		double price_inf = full_mesh(nb_time_steps_, position_inf);
		double price_sup = full_mesh(nb_time_steps_, position_sup);
		price = Mesh::linear_interp_around_spot(log_spots, price_inf, price_sup);
	}

	return price;
}

double Mesh::delta_formula(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots, 
						int const& position_inf, int const& position_sup) const
{
	double dx = Mesh::comp_dx();
	double delta;
	double f_0_sup = full_mesh(nb_time_steps_, position_sup);
	double f_0_inf = full_mesh(nb_time_steps_, position_inf);
	double log_derivative_origin = (f_0_sup - f_0_inf) / (2 * dx);
	int position = position_sup - 1;
	double spot = exp(log_spots(position));
	delta = log_derivative_origin / spot;
	return delta;
}

double Mesh::comp_delta(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const
{
	double dx = Mesh::comp_dx();
	double delta;
	if (nb_spot_steps_ % 2 == 0)
	{
		
		int position_sup = (nb_spot_steps_) / 2 + 1;
		int position_inf = (nb_spot_steps_) / 2 - 1;
		delta = Mesh::delta_formula(full_mesh, log_spots, position_inf, position_sup);

	}
	else
	{
		//linear interpolation between the deltas at log spots around the initial log spot
		int position_sup_sup = (nb_spot_steps_ - 1) / 2 + 2;
		int position_sup_inf = (nb_spot_steps_ - 1) / 2;
		double delta_sup = Mesh::delta_formula(full_mesh, log_spots, position_sup_inf, position_sup_sup);

		int position_inf_sup = (nb_spot_steps_ - 1) / 2 + 1;
		int position_inf_inf = (nb_spot_steps_ - 1) / 2 - 1;
		double delta_inf = Mesh::delta_formula(full_mesh, log_spots, position_inf_inf, position_inf_sup);

		delta = Mesh::linear_interp_around_spot(log_spots, delta_inf, delta_sup);

	}
	return delta;
}

double Mesh::gamma_formula(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots,
	int const& position_inf, int const& position_sup) const
{
	double dx = Mesh::comp_dx();
	double gamma;
	double f_0_sup = full_mesh(nb_time_steps_, position_sup);
	double f_0_inf = full_mesh(nb_time_steps_, position_inf);
	int position = position_sup - 1;
	double price_pos = full_mesh(nb_time_steps_, position);
	double log_derivative_origin = (f_0_sup - f_0_inf) / (2 * dx);
	double second_log_derivative_origin = (f_0_sup - 2 * price_pos + f_0_inf) / pow(dx, 2);
	double spot = exp(log_spots(position));
	gamma = (second_log_derivative_origin - log_derivative_origin) / pow(spot, 2);
	return gamma;
}





double Mesh::comp_gamma(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const
{
	double dx = Mesh::comp_dx();
	double gamma;
	if (nb_spot_steps_ % 2 == 0)
	{
		int position_sup = (nb_spot_steps_) / 2 + 1;
		int position_inf = (nb_spot_steps_) / 2 - 1;
		gamma = Mesh::gamma_formula(full_mesh, log_spots, position_inf, position_sup);
		
	}
	else
	{
		//linear interpolation between the gammas at the two log spots around the level of initial log spot.
		int position_sup_sup = (nb_spot_steps_ - 1) / 2 + 2;
		int position_sup_inf = (nb_spot_steps_ - 1) / 2;
		double gamma_sup = Mesh::gamma_formula(full_mesh, log_spots, position_sup_inf, position_sup_sup);

		int position_inf_sup = (nb_spot_steps_ - 1) / 2 + 1;
		int position_inf_inf = (nb_spot_steps_ - 1) / 2 - 1;
		double gamma_inf = Mesh::gamma_formula(full_mesh, log_spots, position_inf_inf, position_inf_sup);

		gamma = Mesh::linear_interp_around_spot(log_spots, gamma_inf, gamma_sup);
	}
	return gamma;
}


double Mesh::comp_theta_greek(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const
{
	double dt = Mesh::comp_dt();
	double theta;
	if (nb_spot_steps_ % 2 == 0)
	{
		theta = (full_mesh(nb_time_steps_ - 1, nb_spot_steps_ / 2) - full_mesh(nb_time_steps_, nb_spot_steps_ / 2)) / dt;
	}
	else
	{
		//linear interpolation between the two theta at the x_i's around the level of the initial log spot
		int position_inf = (nb_spot_steps_ - 1) / 2;
		int position_sup = (nb_spot_steps_ - 1) / 2 + 1;
		double x_inf = log_spots(position_inf);
		double x_sup = log_spots(position_sup);
		double log_spot = log(initial_spot_);
		double theta_inf = (full_mesh(nb_time_steps_ - 1, (nb_spot_steps_ -1) / 2) - full_mesh(nb_time_steps_, (nb_spot_steps_-1) / 2 )) 
							/ dt;
		double theta_sup = (full_mesh(nb_time_steps_ - 1, (nb_spot_steps_ - 1) / 2 + 1) - full_mesh(nb_time_steps_, (nb_spot_steps_ - 1) / 2 + 1))
						/ dt;
		theta = ((log_spot - x_inf) * theta_sup + (x_sup - log_spot) * theta_inf) / (x_sup - x_inf);
	}
	return theta;
}


double Mesh::comp_vega(const double& original_price, const Eigen::MatrixXd& log_spots) const
{
	bool bumped = true;
	Eigen::MatrixXd new_mesh = Mesh::comp_mesh(log_spots, bumped);
	double new_price = Mesh::extract_price(new_mesh, log_spots);
	double vega = (new_price - original_price) / eps_;
	return vega;
}



