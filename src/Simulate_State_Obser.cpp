#include <RcppArmadillo.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate State and Observations
 //' 
 //' This function simulates the state trajectory and corresponding observations for a given system.
 //' The system dynamics and observation model are described by user-provided functions `f` and `h`.
 //' 
 //' @name Simulate_State_Obser
 //' @alias Simulate_State_Obser
 //' @title Simulate State and Observations
 //' @usage Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed)
 //' @param Dt Numeric. The time step size for the simulation.
 //' @param Ntau Integer. The number of observation points.
 //' @param NtNtau Integer. The total number of time steps in the simulation.
 //' @param f Function. A user-provided function representing the state dynamics.
 //' @param h Function. A user-provided function representing the observation model.
 //' @param Dim Integer. The dimension of the system.
 //' @param seed Nullable Integer. An optional random seed for reproducibility.
 //' @return A list containing:
 //' \item{x}{A matrix of the state trajectory over time.}
 //' \item{y}{A matrix of the observations at specific time points.}
 //' \item{seed}{The random seed used for the simulation.}
 //' @description This function simulates the evolution of a system over time using the provided dynamics and observation models. The state is updated using the dynamics function `f` at each time step, while observations are generated using the function `h` at specific intervals. A random seed can be used to ensure reproducibility.
 //' @examples
 //' Dt <- 0.01
 //' Ntau <- 10
 //' NtNtau <- 100
 //' f <- function(x) { return(x) } # Example dynamics function
 //' h <- function(x) { return(x) } # Example observation function
 //' Dim <- 2
 //' result <- Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed = 123)
 // [[Rcpp::export]]
 Rcpp::List Simulate_State_Obser(double Dt, int Ntau, int NtNtau, Rcpp::Function f, Rcpp::Function h, int Dim, Rcpp::Nullable<int> seed = R_NilValue) {
   arma::mat x(NtNtau, Dim, arma::fill::zeros);
   arma::mat y_Dt(NtNtau, Dim, arma::fill::zeros);
   arma::mat y(Ntau, Dim, arma::fill::zeros);
   double sqrtdT = sqrt(Dt);
   
   // Set random seed for reproducibility
   unsigned int seed_value;
   if (seed.isNull()) {
     seed_value = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count() % 1000000000);
   } else {
     seed_value = static_cast<unsigned int>(Rcpp::as<int>(seed));
   }
   
   // Create random number generator
   std::mt19937 generator(seed_value);
   std::normal_distribution<double> distribution(0.0, 1.0);
   
   Rcpp::Rcout << "Random seed used: " << seed_value << "\n";
   
   // Simulate x
   for (int t = 1; t < NtNtau; ++t) {
     arma::rowvec x_prev = x.row(t - 1);
     Rcpp::NumericVector f_value = Rcpp::as<Rcpp::NumericVector>(f(x_prev));
     arma::rowvec noise(Dim);
     for (int d = 0; d < Dim; ++d) {
       noise[d] = distribution(generator);
     }
     x.row(t) = x_prev + Rcpp::as<arma::rowvec>(f_value) * Dt + sqrtdT * noise;
   }
   
   // Simulate y_Dt
   for (int t = 1; t < NtNtau; ++t) {
     arma::rowvec x_prev = x.row(t);
     Rcpp::NumericVector h_value = Rcpp::as<Rcpp::NumericVector>(h(x_prev));
     arma::rowvec noise(Dim);
     for (int d = 0; d < Dim; ++d) {
       noise[d] = distribution(generator);
     }
     y_Dt.row(t) = y_Dt.row(t - 1) + Rcpp::as<arma::rowvec>(h_value) * Dt + sqrtdT * noise;
   }
   
   // Simulate y
   for (int n = 1; n < Ntau; ++n) {
     int t = n * (NtNtau / Ntau);
     y.row(n) = y_Dt.row(t);
   }
   
   return Rcpp::List::create(Rcpp::Named("x") = x,
                             //Rcpp::Named("y_Dt") = y_Dt,
                             Rcpp::Named("y") = y,
                             Rcpp::Named("seed") = seed_value);
 }
