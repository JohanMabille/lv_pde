#pragma once

#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"
#include "pde_matrix_builder.hpp"
#include "system_solver.hpp"
#include "payoff.hpp"


class Mesh {

public:

	virtual ~Mesh() = default;
	SystemSolver* syst_solv_;
	Payoff* payoff_;

	Mesh(SystemSolver* syst_solv, Payoff* payoff, const int& nb_spot_steps, const int& nb_time_steps, const double& initial_spot,
		const double& K, const double& initial_vol, const double& T, const double& initial_rate, const bool& constant_coeffs, const double& eps);
	const int& get_nb_time_steps() const;
	const int& get_nb_spot_steps() const;
	const double& get_initial_spot() const;

	double comp_dx() const;
	double comp_dt() const;

	Eigen::MatrixXd comp_mesh(const Eigen::MatrixXd& log_spots, const bool bumped=false) const;

	Eigen::MatrixXd mesh_maturity() const;
	Eigen::MatrixXd log_spot_prices() const;

	double extract_price(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const;

	
	double comp_delta(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const;
	double comp_gamma(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const;
	double comp_theta_greek(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots) const;
	double comp_vega(const double& original_price, const Eigen::MatrixXd& log_spots) const;

	double get_K() const;
	double get_vol() const;
	double get_T() const;
	double get_rate() const;


private:

	bool constant_coeffs_;
	double K_;
	double initial_vol_;
	double T_;
	double initial_rate_;
	int nb_time_steps_;
	int nb_spot_steps_;
	double initial_spot_;
	double eps_;

	double linear_interp_around_spot(const Eigen::MatrixXd& log_spots, const double& y_inf, const double& y_sup) const;
	double delta_formula(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots, int const& position_inf, int const& position_sup) const;
	double gamma_formula(const Eigen::MatrixXd& full_mesh, const Eigen::MatrixXd& log_spots, int const& position_inf, int const& position_sup) const;
	
};

