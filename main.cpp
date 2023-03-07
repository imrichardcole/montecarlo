#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

typedef struct Trade {

  Trade(int id, double option_price, double strike_price, double risk_free_rate, double vol, double time_to_expiry) : 
    id(id),
    option_price(option_price), 
    strike_price(strike_price),
    risk_free_rate(risk_free_rate),
    vol(vol),
    time_to_expiry(time_to_expiry) {};

  const int id;
  const double option_price;
  const double strike_price;
  const double risk_free_rate;
  const double vol;
  const double time_to_expiry;

} Trade;

double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}

double black_scholes_call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return S * norm_cdf(d_j(1, S, K, r, v, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, v, T));
}

double black_scholes_put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return -S * norm_cdf(-d_j(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);
  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

double monte_carlo_call_price(const int& simulations, Trade &trade) {
  double S_adjust = trade.strike_price * exp(trade.time_to_expiry * ( trade.risk_free_rate - 0.5 * trade.vol * trade.vol));
  double S_cur = 0.0;
  double payoff_sum = 0.0;

  for (int i=0; i < simulations; i++) {
    double gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * exp(sqrt(trade.vol * trade.vol * trade.time_to_expiry) * gauss_bm);
    payoff_sum += std::max(S_cur - trade.strike_price, 0.0);
  }

  return (payoff_sum / static_cast<double>(simulations)) * exp(-trade.risk_free_rate * trade.time_to_expiry);
}

double monte_carlo_put_price(const int& simulations, Trade& trade) {
  double S_adjust = trade.strike_price * exp(trade.time_to_expiry * ( trade.risk_free_rate - 0.5 * trade.vol * trade.vol));
  double S_cur = 0.0;
  double payoff_sum = 0.0;

  for (int i=0; i < simulations; i++) {
    double gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * exp(sqrt(trade.vol * trade.vol * trade.time_to_expiry) * gauss_bm);
    payoff_sum += std::max(trade.strike_price - S_cur, 0.0);
  }

  return (payoff_sum / static_cast<double>(simulations)) * exp(-trade.risk_free_rate * trade.time_to_expiry);
}


int main(int argc, char const *argv[])
{
  const std::vector<Trade> trades = {
    {1, 100, 100, 0.05, 0.2, 1},
    {2, 100, 101, 0.04, 0.2, 2},
    {3, 120, 99, 0.06, 0.2, 5},
    {4, 14, 10, 0.045, 0.2, 10}
  };

  const int simulation_count = 1000000;

  for(auto trade : trades) {
    const double put_price = monte_carlo_put_price(simulation_count, trade);
    const double bs_put_price = black_scholes_put_price(trade.option_price, trade.strike_price, trade.risk_free_rate, trade.vol, trade.time_to_expiry);

    const double call_price = monte_carlo_call_price(simulation_count, trade);
    const double bs_call_price = black_scholes_call_price(trade.option_price, trade.strike_price, trade.risk_free_rate, trade.vol, trade.time_to_expiry);
  
    std::cout << "******** Trade " << trade.id << " ********" << std::endl;
    std::cout << " MC put price = " << put_price << std::endl;
    std::cout << " BS put price = " << bs_put_price << std::endl;
    std::cout << " MC call price = " << call_price << std::endl;
    std::cout << " BS call price = " << bs_call_price << std::endl;
    std::cout << "************************" << std::endl;
  }

  return 0;
}
