#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP

#include <vector>


    double vanilla_payoff(double fwd, double strike, bool is_call);
    double bs_time_value(double fwd, double strike, double volatility, double maturity);
    double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

    double delta_call(const double spot, const double K, const double rate, const double vol, const double T);
    double gamma_call(const double spot, const double K, const double rate, const double vol, const double T);
    double vega_call(const double spot, const double K, const double rate, const double vol, const double T);
    double theta_call(const double spot, const double K, const double rate, const double vol, const double T);

    double delta_put(const double spot, const double K, const double rate, const double vol, const double T);
    double gamma_put(const double spot, const double K, const double rate, const double vol, const double T);
    double vega_put(const double spot, const double K, const double rate, const double vol, const double T);
    double theta_put(const double spot, const double K, const double rate, const double vol, const double T);

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);


#endif
