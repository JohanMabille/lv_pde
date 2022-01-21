#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>
#include <execution>
#include "payoff.hpp"
#include "Option.hpp"
#include "Financial_PDE.hpp"
#include "boundary_conditions.hpp"
#include "pde_matrix_builder.hpp"
#include <Eigen/Dense>
#include "system_solver.hpp"
#include "mesh.hpp"
#include "closed_form.hpp"


void test_boundaryconditions()
{

}

int main() {
	double K = 10.0;
	double r = 0.05;
	double vol = 0.1;
	double T = 1.0;

	//double dx = 0.1;
	//double dt = 0.1;
	double theta = 0.5;
	
	int nb_time_steps = 252;
	int nb_spot_steps = 90;
	double initial_spot = 8.72;

	Payoff* call_payoff = new PutPayoff(K);
	Option* call = new VanillaOption(call_payoff, K, vol, T, r);
	Financial_PDE* bs_pde = new BS_PDE(call);

	BoundaryConditions* bound_x0 = new Boundaryx0(bs_pde, theta);
	BoundaryConditions* bound_xN = new BoundaryxN(bs_pde, theta);

	// Problème pour créer un objet de type PDEMatrixBuilder... Tout le reste au-dessus marche bien...

	PDEMatrixBuilder* mat_build = new PDEMatrixBuilder(bound_x0, bound_xN);

	// TEST FOR MATRIX COMPUTATION

	/*
	Eigen::MatrixXd k = mat_build->comp_k();
	std::cout << k << std::endl;

	Eigen::MatrixXd ktilde = mat_build->comp_ktilde();
	std::cout << ktilde << std::endl;

	Eigen::MatrixXd fn1 = Eigen::MatrixXd::Random(nb_time_steps + 1, 1);
	Eigen::MatrixXd fn = syst_solv->get_fn(fn1);

	std::cout << fn1 << std::endl;
	std::cout << fn << std::endl;

	*/

	SystemSolver* syst_solv = new SystemSolver(mat_build);
	Mesh* mesh = new Mesh(syst_solv, nb_spot_steps, nb_time_steps, initial_spot);

	Eigen::MatrixXd log_spots = mesh->log_spot_prices();
	Eigen::MatrixXd mesh_prices = mesh->comp_mesh();
	//std::cout << mesh_prices << std::endl;
	double price = mesh->extract_price(mesh_prices, log_spots);
	std::cout << mesh_prices << std::endl;
	std::cout << price << std::endl;
	/*
	double delta = mesh->comp_delta(mesh_prices, log_spots);
	std::cout << delta << std::endl;
	double gamma = mesh->comp_gamma(mesh_prices, log_spots, price);
	std::cout << gamma << std::endl;
	*/


	// RAJOUTER LES DELETE POINTEURS!
	delete mesh;
	mesh = NULL;
	delete mat_build;
	mat_build = NULL;
	delete bound_xN;
	bound_xN = NULL;
	delete bound_x0;
	bound_x0 = NULL;
	delete bs_pde;
	bs_pde = NULL;
	delete call;
	call = NULL;
	delete call_payoff;
	call_payoff = NULL;





	//double closed_form_price = bs_price()






	//Eigen::MatrixXd result_mesh = mesh->comp_mesh();
	//Eigen::MatrixXd res = Eigen::MatrixXd::Zero(nb_time_steps + 1, nb_spot_steps + 1);

	//res.row(0) = mesh->mesh_maturity();
	//Eigen::MatrixXd fn1 = Eigen::MatrixXd::Zero(nb_spot_steps +1, 1);


	/*
	Eigen::MatrixXd log_spots = mesh->log_spot_prices();
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1, nb_spot_steps + 1);
	for (int i = 0; i <= nb_spot_steps; i++)
	{
		double log_spot = log_spots(0, i);
		std::cout << log_spot << std::endl;
		res(0, i) = syst_solv->mat_build_->boundaries_x0_->pde_->Opt_->m_payoff->operator()(exp(log_spot));
		std::cout << res(0,i) << std::endl;
	}

	*/




	return 0;
}