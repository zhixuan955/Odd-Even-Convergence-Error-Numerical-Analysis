# Odd-Even-Convergence-Error-Numerical-Analysis
A numerical investigation and statistical analysis of the Odd/Even effect in the Cox-Ross-Rubinstein binomial tree as compared against the closed form Black Scholes Model. 

Options pricing through a Cox-Ross-Rubinstein binomial tree involves partitioning the option lifetime into discrete steps dt = T/N (where T is time at maturity and N is the total number of steps). It is intuitive to conclude that increasing N (therefore reducing the size of each time increment) constitutes a smooth and consistent convergence of the price to that modelled by the Black-Scholes model. 

However, this convergence is instead not smooth and involves oscillating up and down depending on whether N is an odd or even number. My goal with this investigation is to show that these oscillations can be mapped predictably onto a function of where the strike sits among the terminal nodes. Furthermore, this script aims to show that the error can be reduced through subtracting a predicted oscillating term at no extra computational cost.

# Headline results
All results are for a European call with:
Initial stock price ($$S0$$) = 100
Time in years ($$T$$) = 1
Risk-free rate ($$r$$) = 5%
Volatility ($$sigma$$) = 20%
Variable strike price ($$K$$) = 70 to 130
Variable number of steps ($$N$$) = 1000 - 1999

1. The oscillation is a quadratic in the fractional node position (xi) of the strike’s node at maturity (t = T), with the equation:
$$ N x (tree - BS) = a(K) + b(K) x B2(xi) $$

where $$ B2(xi) = xi^2 - xi + 1/6 $$ and B2(x) $$ is the second Bernoulli equation
Iterating through 60 strike prices (excluding $$ K = S0 $$ ) and applying linear least squares regression produced a coefficient of determination ( $$ R^2 #$$ ) > 0.9989 with median 0.9999

3. 
