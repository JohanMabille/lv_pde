#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>
#include <execution>


#include "payoff.hpp"
#include "Financial_PDE.hpp"
#include "boundary_conditions.hpp"
#include "pde_matrix_builder.hpp"
#include "eigen-3.4.0/Eigen/Dense"
#include "system_solver.hpp"
#include "mesh.hpp"
#include "closed_form.hpp"
#include "volatility_diffusion.hpp"
#include "rate_diffusion.hpp"
#include "closed_form.hpp"


void test_pricer_greeks(double K, double initial_r, double initial_vol, double T, double theta, double nb_time_steps, double nb_spot_steps,
	double initial_spot, double constant_coeffs, double eps, Payoff* payoff, VolatilityDiffusion* vol_diff, RateDiffusion* rate,
	Financial_PDE* bs_pde, BoundaryConditions* bound_x0, BoundaryConditions* bound_xN, PDEMatrixBuilder* mat_build, SystemSolver* syst_solv,
	Mesh* mesh, int option_type)
{
	Eigen::MatrixXd log_spots = mesh->log_spot_prices();
	Eigen::MatrixXd mesh_prices = mesh->comp_mesh(log_spots);

	// Mesh values
	double price = mesh->extract_price(mesh_prices, log_spots);
	std::cout << "Mesh price: " << price << std::endl;
	double delta = mesh->comp_delta(mesh_prices, log_spots);
	std::cout << "Mesh Delta: " << delta << std::endl;
	double gamma = mesh->comp_gamma(mesh_prices, log_spots);
	std::cout << "Mesh Gamma: " << gamma << std::endl;
	double theta_greek = mesh->comp_theta_greek(mesh_prices, log_spots);
	std::cout << "Mesh Theta: " << theta_greek << std::endl;
	double vega = mesh->comp_vega(price, log_spots);
	std::cout << "Mesh Vega: " << vega << std::endl;

        // If you provide the possibility to price a call or a put, you should handle it
        // here too.
        double price_explicit, delta_explicit, gamma_explicit, theta_explicit, vega_explicit;
	// Closed form values
	// CALL Option:
	if (option_type == 1)
        {
            price_explicit = bs_price(initial_spot * exp(initial_r * T), K, initial_vol, T, true);
            std::cout << "closed form price: " << price_explicit << std::endl;
            delta_explicit = delta_call(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form delta: " << delta_explicit << std::endl;
            gamma_explicit = gamma_call(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form gamma: " << gamma_explicit << std::endl;
            theta_explicit = theta_call(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form theta: " << theta_explicit << std::endl;
            vega_explicit = vega_call(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form vega: " << vega_explicit << std::endl;
        }
	// PUT Option:
        else
        {
            price_explicit = bs_price(initial_spot * exp(initial_r * T), K, initial_vol, T, false);
            std::cout << "closed form price: " << price_explicit << std::endl;
            delta_explicit = delta_put(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form delta: " << delta_explicit << std::endl;
            gamma_explicit = gamma_put(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form gamma: " << gamma_explicit << std::endl;
            theta_explicit = theta_put(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form theta: " << theta_explicit << std::endl;
            vega_explicit = vega_put(initial_spot, K, initial_r, initial_vol, T);
            std::cout << "closed form vega: " << vega_explicit << std::endl;
        }
	
	

	// comparing the results with the two methods (error expressed in percentage):

	double price_error = abs((price - price_explicit) / price_explicit);
	std::cout << "price percentage error: " << price_error << std::endl;
	double delta_error = abs((delta - delta_explicit) / delta_explicit);
	std::cout << "delta percentage error: " << delta_error << std::endl;
	double gamma_error = abs((gamma - gamma_explicit) / gamma_explicit);
	std::cout << "gamma percentage error: " << gamma_error << std::endl;
	double theta_error = abs((theta_greek - theta_explicit) / theta_explicit);
	std::cout << "theta percentage error: " << theta_error << std::endl;
	double vega_error = abs((vega - vega_explicit) / vega_explicit);
	std::cout << "vega percentage error: " << vega_error << std::endl;
}

int main() {
	
	//Option Parameters:
	/*std::cout << "Option parameters" << std::endl;
	double K;
	std::cout << "please enter the strike" << std::endl;
	std::cin >> K;
	double initial_r;
	std::cout << "please enter the initial interest rate" << std::endl;
	std::cin >> initial_r;
	double initial_vol;
	std::cout << "please enter the initial volatility" << std::endl;
	std::cin >> initial_vol;
	double T;
	std::cout << "please enter the maturity in years" << std::endl;
	std::cin >> T;

	//Theta scheme parameter:
	double theta;
	std::cout << "please enter the parameter of the Theta scheme (between zero and one)" << std::endl;
	std::cin >> theta;

	// Mesh parameters:
	std::cout << "Mesh parameters" << std::endl;
	int nb_time_steps;
	std::cout << "please enter the number of time steps in the mesh" << std::endl;
	std::cin >> nb_time_steps;
	int nb_spot_steps;
	std::cout << "please enter the number of spot steps in the mesh" << std::endl;
	std::cin >> nb_spot_steps;
	double initial_spot;
	std::cout << "please enter the spot price" << std::endl;
	std::cin >> initial_spot;

	int option_type;
	std::cout << "Enter 0 if the option is a Put, enter 1 if it's a call" << std::endl;
	std::cin >> option_type;

        // This is WRONG: you don't capture boolean values with strings that you store in double
        // The result will never be what you expect
	double constant_coeffs;
	std::cout << "if the model uses constant diffusion coefficients enter true, else enter false" << std::endl;
	std::cin >>constant_coeffs;*/

        double K = 150.;
        double initial_spot = 100.;
        double initial_r = 0.0;
        double initial_vol = 0.2;
        double T = 1.;
        double theta = 0.5;
        int nb_time_steps = 365;
        int nb_spot_steps = 501;
        int option_type = 1;
        // pricing with constant_coeffs = false never completes...
        bool constant_coeffs = true;

        std::cout << "constant_coeffs = " << constant_coeffs << std::endl;
	double eps = 0.001;

	
	Payoff* payoff;
	if (option_type == 1)
	{
		payoff = new CallPayoff(K);
	}
	else
	{
		payoff = new PutPayoff(K);
	}
	
        // Dynamic allocation is not required here, you could allocate on the stack
        // and pass the objects by reference to test_pricer_greeks (so that polymorphism
        // works as expected)
	VolatilityDiffusion* vol_diff = new ConstantVolatility();
	RateDiffusion* rate_diff = new ConstantRate();
	Financial_PDE* bs_pde = new BS_PDE(vol_diff, rate_diff);
	BoundaryConditions* bound_x0 = new Boundaryx0(bs_pde, theta);
	BoundaryConditions* bound_xN = new BoundaryxN(bs_pde, theta);
	PDEMatrixBuilder* mat_build = new PDEMatrixBuilder(bound_x0, bound_xN);
	SystemSolver* syst_solv = new SystemSolver(mat_build);
	Mesh* mesh = new Mesh(syst_solv, payoff, nb_spot_steps, nb_time_steps, initial_spot, K, initial_vol, T, initial_r, constant_coeffs, eps);

	test_pricer_greeks(K, initial_r, initial_vol, T, theta, nb_time_steps, nb_spot_steps, initial_spot, constant_coeffs, eps, payoff, vol_diff,
		rate_diff, bs_pde, bound_x0, bound_xN, mat_build, syst_solv, mesh, option_type);

	delete mesh;
	mesh = nullptr;
        // Where do you delete syst_solv?
	delete mat_build;
	mat_build = nullptr;
        // The following leads to a crash
	delete bound_xN;
	bound_xN = nullptr;
	delete bound_x0;
	bound_x0 = nullptr;
	delete bs_pde;
	bs_pde = nullptr;
	delete payoff;
	payoff = nullptr;


	return 0;
}
