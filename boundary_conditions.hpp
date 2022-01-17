#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "Financial_PDE.hpp"

class BoundaryConditions
{
public:
	virtual ~BoundaryConditions() = default;
	//BoundaryConditions();

	Financial_PDE* pde_; // had to pass the financial_PDE pointer to public for it to be accessible from a pointer to BoundaryCondition object.

	const double& get_theta() const;
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt) const = 0;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt) const = 0;

protected:

	
	BoundaryConditions(Financial_PDE* pde, const double& theta);


private:

	double theta_;

};

class Boundaryx0: public BoundaryConditions 
{

public:

	Boundaryx0(Financial_PDE* pde, const double& theta);
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt) const override;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt) const override;

protected:

	
	

private:

	/*
	double gamma_x0(double& dx) const;
	double vega_x0(double& dx) const;
	double mu_x0(double& dx) const;
	*/
	double gamma_x0(const double& dx, const double& dt) const;
	double vega_x0(const double& dx, const double& dt) const;
	double mu_x0(const double& dx, const double& dt) const;
};

class BoundaryxN : public BoundaryConditions
{

public:

	BoundaryxN(Financial_PDE* pde, const double& theta);
	virtual Eigen::MatrixXd coeff_fn(const double& dx, const double& dt) const override;
	virtual Eigen::MatrixXd coeff_fn1(const double& dx, const double& dt) const override;

protected:

	

private:

	/*
	double gamma_x0(double& dx) const;
	double vega_x0(double& dx) const;
	double mu_x0(double& dx) const;
	*/
	double gamma_xN(const double& dx, const double& dt) const;
	double vega_xN(const double& dx, const double& dt) const;
	double mu_xN(const double& dx, const double& dt) const;
};


/*
3 classes filles � faire:
	- Boundaryx0: gauche de la mesh quand log(s)=x_0
	- BoundaryxN: droite de la mesh quand log(s)=x_N
	- BoundaryT: haut de la mesh quand t=T (payoff � maturit� selon exo(log(s)))
qui correspondent aux trois boundaries connues.

Il faudrait que �a h�rite d'un pointeur de type option encore, pour avoir les donn�es de r, sigma, T, et du payoff.

EN FAIT, juste besoin d'h�riter d'un pointeur de type FinancialPDE. Puis, on pourra r�cup�rer l'option, le payoff etc et la valeur
des coefficients au temps t. Mettre un param�tre t pour les fonctions sur les coefficients peut-�tre je pense dans FinancialPDE, dans un 
souci de g�n�ralit�.

Les param�tres n�cessaires autre que ce qui est d�j� contenu dans le pointeur vers un objet de type Financial_PDE sont simplement dx et dt.
Il faut cr�er des fonctions qui ressortent des vecteurs de dimension trois en vrai avec vect=(A_0, B_0, C_0) par exemple pour x0.
On va avoir un vecteur qui correspond � f^n et un � f^(n+1).
Donc, seulement besoin de deux fonctions en vrai.

ATTENTION, faire plus de fonctions en fait. Il faudrait faire des sous fonctions qui calculent ce que j'ai not� gamma, vega, et mu sur la 
feuille.

ATTENTION, il y a en plus besoin du param�tre Theta qui doit pouvoir �tre modifi� par l'utilisateur je pense.
*/