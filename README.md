# Odd-Even-Convergence-Error-Numerical-Analysis
A numerical investigation and statistical analysis of the Odd/Even effect in the Cox-Ross-Rubinstein binomial tree as compared against the closed form Black Scholes Model. 

![Scatter plot of scaled error against \xi, scatter plot of scaled error divided by amplitude against $$B_2(\xi)$$](OE_figure.png)

Error convergence in a Cox-Ross-Rubinstein binomial tree is not smooth and involves oscillations up and down depending on whether N is an odd or even number. My goal with this investigation is to show that these oscillations can be mapped predictably onto a function of where the strike sits among the terminal nodes. Furthermore, this script aims to show that the error can be reduced through subtracting a predicted oscillating term at no extra computational cost.

# Headline results
All results are for a European call with:
Initial stock price ($$S0$$) = 100
Time in years ($$T$$) = 1
Risk-free rate ($$r$$) = 5%
Volatility ($$sigma$$) = 20%
Variable strike price ($$K$$) = 70 to 130
Variable number of steps ($$N$$) = 1000 - 1999

- The oscillation is a quadratic in the fractional node position (xi) of the strike’s node at maturity (t = T), with the equation: $$N \times (\text{tree} - \text{BS}) = a(K) + b(K) \times B_2(\xi) \qquad$$ where $$B_2(\xi) = (\xi)^2 - \xi + 1/6$$ and $$B2(x)$$ is the second Bernoulli equation. Iterating through 60 strike prices (excluding $$ K = S0 $$ ) and applying linear least squares regression produced a coefficient of determination ( $$ R^2 #$$ ) > 0.9989 with median 0.9999

- The amplitude matches a formula from a heuristic argument. The fitted $$b(K)$$ equals $$-2\sigma \times \sqrt{T} \times g(K)$$ where $$g(K) = K \exp{-r \times T} \phi(d_2)$$, to within 0.1%. (b/g ratio between -0.4003 and -0.3998 and theoretical ratio equal to -0.4)

- Subtracting the predicted term lowers the error at no extra computational cost. This led to median absolute error to shrink 2.25x and the maximum absolute error to shrink 3.59x (table below)

# Background
A call option is a derivative contract that derives its value from the performance of an underlying asset. For a premium (upfront price) to be determined, a call option gives the buyer the right, but not the obligation to buy the asset at a predetermined price (the strike price) on or before a the expiration date. The Black Scholes model is a financial model used to estimate the theoretical exact value of financial options. For the purpose of benchmarking the Cox-Ross-Rubinstein against the closed form Black Scholes, this investigation will focus on European Call Options, which allow the buyer to sell only on the expiration date. 

Options pricing through a Cox-Ross-Rubinstein binomial tree involves partitioning the option lifetime into discrete steps $$dt = T/N$$ (where $$T$$ is time at maturity and N is the total number of steps). After each time step, the stock price can either go up or go down an equivalent amounts. At the expiration date, the final price of the stock will be the expected payoff of the option if it were to be exercised by the buyer. The premium to buy this call option follows the formula:
$$C = discounted expected payoff under risk-neutral probabilities$$ 

It is intuitive to conclude that increasing N (therefore reducing the size of each time increment) constitutes a smooth and consistent convergence of the price to that modelled by the Black-Scholes model. However, this does not happen due to the odd/even effect, where the error between option price computed by the tree and that of Black-Scholes does not shrink smoothly as N grows, but instead alternates depending on whether N is odd or even.
  
# Methodology
**Models** A Cox-Ross_Rubinstein binomial tree where the up factor ($$u = exp{sigma \sqrt{T}}) and down factor ($$d = 1/u$$) allow the tree to recombine. The model is implemented in the script to adopt two different methods of computing the option price. First is backward induction through the tree, second is a faster equivalent that sums the terminal payoff weighted by the binomial probabilities of $$j$$ number of up movements that occurred during the lifetime of the option. Validation check within the script returned a difference of 2.1e-12 between the two pricing methods across 9 test cases (K = 70, 100, 130; N = 101, 500, 1000). In every practical sense, they agree completely. The script moves forward with the fast version. The Black-Scholes model is used to compute the exact option price that will be compared against in subsequent tests. 
