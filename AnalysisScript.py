#-----------------------------------------------------------------------------------------------------------------------
# Import modules
#-----------------------------------------------------------------------------------------------------------------------
import numpy as np
import scipy.stats as stats
import csv
import matplotlib.pyplot as plt
import pandas as pd


#-----------------------------------------------------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------------------------------------------------
S0 = 100        # initial stock price
T = 1           # expiration time (in years)
r = 0.05        # risk-free interest rate
sigma = 0.2     # volatility

strike_prices = np.arange(70, 130.01, 1.0)  # K = 70, 71, ... , 130
N_values = np.arange(1000, 2000)                 # N = 1000, 1001, ... , 1999


#-----------------------------------------------------------------------------------------------------------------------
# Black-Scholes Model (closed form)
#-----------------------------------------------------------------------------------------------------------------------
def black_scholes(S0, K, T, r, sigma, opttype="Call"):

    # computing d1 and d2
    d1 = (np.log(S0 / K) + (r + (sigma**2 / 2)) * T)/(sigma * np.sqrt(T))
    d2 = d1 - (sigma * np.sqrt(T))
    discount = np.exp(-r * T)

    if opttype == "Call":
        # initiate call option
        return S0 * stats.norm.cdf(d1) - K * discount * stats.norm.cdf(d2)
    else:
        # initiate put option
        return K * discount * stats.norm.cdf(-d2) - S0 * stats.norm.cdf(-d1)


#-----------------------------------------------------------------------------------------------------------------------
# Baseline Binomial Tree Model using backwards induction
#-----------------------------------------------------------------------------------------------------------------------
def european_binomial_tree(S0, K, T, N, r, sigma, opttype="Call"):

    dt = T / N  # time increment size
    discount = np.exp(-r * dt)

    # computing up and down factors
    u = np.exp(sigma * (dt ** 0.5))                         # d = 1/u

    risk_neutral_prob = ((1 / discount) - u ** (-1)) / (u - u ** (-1))
    # computing stock price with the number of up and down movements
    stock = S0 * (u ** (-1)) ** (np.arange(N, -1, -1)) * u ** (np.arange(0, N + 1, 1))

    if opttype == 'Call':
        option = np.maximum(stock - K, np.zeros(N + 1))
    else:
        option = np.maximum(K - stock, np.zeros(N + 1))

    # backwards induction and returns present value of option price
    for i in np.arange(N, 0, -1):
        option = discount * (risk_neutral_prob * option[1:i + 1] + (1 - risk_neutral_prob) * option[0:i])

    return option[0]


#-----------------------------------------------------------------------------------------------------------------------
# Fast Binomial Tree Model
#-----------------------------------------------------------------------------------------------------------------------
def european_binomial_tree_fast(S0, K, T, N, r, sigma, opttype="Call"):

    dt = T / N  # time increment size
    discount = np.exp(-r * dt)

    # computing up and down factors
    u = np.exp(sigma * (dt ** 0.5))                          # d = 1/u

    risk_neutral_prob = ((1 / discount) - u**(-1)) / (u - u**(-1))
    j = np.arange(N + 1)                                     # number of up moves in total
    stock = S0 * np.exp((2 * j - N) * sigma * np.sqrt(dt))   # S0 * u**j * d**(N-j)

    if opttype == "Call":
        payoff = np.maximum(stock - K, 0)
    else:
        payoff = np.maximum(K - stock, 0)

    # Computes backward induction result in one go
    return np.exp(-r * T) * np.sum(stats.binom.pmf(j, N, risk_neutral_prob) * payoff)
        # stats.binom.pmf(j, N, risk_neutral_prob) - probability of exactly j up movements (j defined as np array)

#-----------------------------------------------------------------------------------------------------------------------
# Validating fast model against baseline
#-----------------------------------------------------------------------------------------------------------------------
print("Validating fast model against baseline model...")
max_diff = 0.0

# testing with one OTM, one ATM and one ITM scenario
for K in (70, 100, 130):
    for N in (101, 500, 1000):
        diff = abs(european_binomial_tree(S0, K, T, N, r, sigma, opttype="Call")
                   - european_binomial_tree_fast(S0, K, T, N, r, sigma, opttype="Call")
                   )
        max_diff = max(max_diff, diff)

# returns error if error too great -> used for debugging only
assert max_diff < 1e-9, f"Validation failed: max difference {max_diff:.2e}"
print(f"Validation passed: maximum difference across all checks is {max_diff:.2e}")


#-----------------------------------------------------------------------------------------------------------------------
# Set up for hypothesis testing
#-----------------------------------------------------------------------------------------------------------------------
def node_fraction(S0, K, T, N, sigma):
    dt = T/N
    j_star = 0.5 * (N + (np.log(K/S0) / (sigma * np.sqrt(dt))))
    # this calculates fractional node position
    return np.round(j_star, 8) % 1.0                      # rounds j_star first, then takes the fractional part
        # If node lands exactly on K, then xi = 0
        # If node straddles K, then xi != 0

def scaled_error(S0, K, T, N, r, sigma):
    error = (european_binomial_tree_fast(S0, K, T, N, r, sigma, opttype="Call")
             - black_scholes(S0, K, T, r, sigma, opttype="Call"))
    return error * N    # scale to N

def second_bernoulli(x):
    # Testing if this parabola in xi fits the curve
    return x**2 - x + 1/6

def strike_price_payment(S0, K, T, r, sigma):
    # d1 - sigma*sqrt(T) = d2 = (ln(S0/K) + (r - sigma^2/2)T) / (sigma*sqrt(T))
    d2 = (np.log(S0 / K) + (r + (sigma**2 / 2)) * T) / (sigma * np.sqrt(T)) - sigma * np.sqrt(T)
    discount = np.exp(-r * T)
    return K * discount * stats.norm.pdf(d2)

# ----------------------------------------------------------------------------------------------------------------------
# Correction method setup
# ----------------------------------------------------------------------------------------------------------------------
def oscillation_term(S0, K, T, N, r, sigma):
    return (-2 * sigma * np.sqrt(T) * strike_price_payment(S0, K, T, r, sigma)
            * second_bernoulli(node_fraction(S0, K, T, N, sigma)) / N)

# the main method of correction
def corrected_price(S0, K, T, N, r, sigma):
    return european_binomial_tree_fast(S0, K, T, N, r, sigma) - oscillation_term(S0, K, T, N, r, sigma)


# ----------------------------------------------------------------------------------------------------------------------
# Testing hypothesis 1: Error depends on xi
# ----------------------------------------------------------------------------------------------------------------------
relation_results = []
for K in strike_prices:
    for N in N_values:
        xi = node_fraction(S0, K, T, N, sigma)

        # appending results to dictionary
        relation_results.append({'K': K,
                                 'N':N,
                                 "xi": xi,
                                 "scaled_error": scaled_error(S0, K, T, N, r, sigma)
                                 })

# writing the .csv file with data from the dictionary
with open("OE_relation_results.csv", 'w', newline='') as csvfile:
    w = csv.DictWriter(csvfile, ['K', 'N', 'xi', 'scaled_error'])
    w.writeheader()
    w.writerows(relation_results)

# used for data presentation later
xi_data = np.array([x["xi"] for x in relation_results])
error_data = np.array([x["scaled_error"] for x in relation_results])

# completion feedback --------------------------------------------------------------------------------------------------
print("\nComputation process complete, written relation data between scaled error and xi into results.csv")
print(f"  Scaled error: min = {error_data.min():.2e}, max = {error_data.max():.2e}")
print(f"  ξ: min = {xi_data.min():.2e}, max = {xi_data.max():.2e}")


# ----------------------------------------------------------------------------------------------------------------------
# Testing hypothesis 2: Scaled error follows curve a(K) + b(K) * B2(xi)
# ----------------------------------------------------------------------------------------------------------------------
result_dataframe = pd.DataFrame(relation_results)             # creating data frame from previous data
                                                              # for subsequent testing
fit_results = []

for K, group in result_dataframe.groupby("K"):    # grouping rows by each unique value of K
    if K == S0:
        continue

    # setting up linear least squares regression
    matrix = np.c_[np.ones(len(group)), second_bernoulli(group['xi'].to_numpy())]
    (a, b), *_ = np.linalg.lstsq(matrix, group['scaled_error'].to_numpy(), rcond=None)
    err = group['scaled_error'].to_numpy()

    # computing R^2 (coefficient_of_determination) for goodness of fit
    ss_res = float(np.sum( (err - matrix @ [a, b])**2 ))
    ss_tot = float(np.sum( (err - err.mean())**2 ))
    coefficient_of_determination = 1 - ss_res/ss_tot

    # appending results to dictionary
    fit_results.append({"K": K,
                        "a": a,
                        "b": b,
                        "b/g_ratio": b / strike_price_payment(S0, K, T, r, sigma),
                        "R^2": coefficient_of_determination
                        })

# writing to the .csv file
with open("OE_fit_results.csv", 'w', newline='') as csvfile:
    w = csv.DictWriter(csvfile, ['K', 'a', 'b', 'b/g_ratio', 'R^2'])
    w.writeheader()
    w.writerows(fit_results)

r2 = np.array([x["R^2"] for x in fit_results])
ratios = np.array([x["b/g_ratio"] for x in fit_results])

# completion feedback --------------------------------------------------------------------------------------------------
print("\nComputation process complete, written result of linear least squares regression into results.csv")

print(f"  strikes fitted: {(len(fit_results))} (K = {(strike_prices[0])} to {(strike_prices[-1])}, K=S0 excluded)")
print(f"  R^2 (goodness of fit) of the quadratic-in-xi fit : min = {r2.min():.4f} median = {np.median(r2):.4f}")

print(f"  b/g: min {ratios.min():.4f}   "         # minimum ratio value
      f"  median {np.median(ratios):.4f}   "      # median ratio value
      f"  max {ratios.max():.4f}   "              # maximum ratio value
      f"  (theory {-2 * sigma * np.sqrt(T)})"     # theoretical value
      )
    # these values should closely match the theoretical value of -0.4 -> validate hypothesis

# investigating at the money scenario
atm_scenario = result_dataframe[result_dataframe["K"] == S0]

even_atm_err = atm_scenario[atm_scenario["N"] % 2 == 0]["scaled_error"].mean()
odd_atm_err = atm_scenario[atm_scenario["N"] % 2 == 1]["scaled_error"].mean()

print(f"At K = S0: even N -> N*error = {even_atm_err:.4f}   odd N -> {odd_atm_err:.4f}")

# ----------------------------------------------------------------------------------------------------------------------
# Evaluating correction method
# ----------------------------------------------------------------------------------------------------------------------
tree_array = np.zeros((len(strike_prices), len(N_values)))
corrected_array = np.zeros((len(strike_prices), len(N_values)))

for a_, K in enumerate(strike_prices):
    # computing black scholes price beforehand to save on computational power
    bs_result = black_scholes(S0, K, T, r, sigma, opttype="Call")

    for b_, N in enumerate(N_values):
        # computing values for normal tree and fast tree
        tree_result = european_binomial_tree_fast(S0, K, T, N, r, sigma, opttype="Call")
        corrected_result = corrected_price(S0, K, T, N, r, sigma)

        # formulating the matrices
        tree_array[a_, b_] = abs(tree_result - bs_result)
        corrected_array[a_, b_] = abs(corrected_result - bs_result)

# Stats on median, 95th percentile and maximum error
tree_stats = (np.median(tree_array), np.percentile(tree_array, 95), np.max(tree_array))
corrected_stats = (np.median(corrected_array), np.percentile(corrected_array, 95), np.max(corrected_array))

# Quantifying the improvement over the base binomial tree model:
result = []
for t, c in zip(tree_stats, corrected_stats):
    ratio = t / c
    result.append(ratio)
error_shrinkage = tuple(result)

# completion feedback --------------------------------------------------------------------------------------------------
print("\nComputation process complete, evaluated correction method.")
print(f"  {'Model':<14} {'Median':>10} {'P95':>10} {'Max':>10}")
print(f"  {'Plain Tree':<14} {tree_stats[0]:>10.2e} {tree_stats[1]:>10.2e} {tree_stats[2]:>10.2e}")
print(f"  {'Corrected':<14} {corrected_stats[0]:>10.2e} {corrected_stats[1]:>10.2e} {corrected_stats[2]:>10.2e}")
print(f"  {'Gain':<14} {error_shrinkage[0]:>9.2f}x {error_shrinkage[1]:>9.2f}x {error_shrinkage[2]:>9.2f}x")


#-----------------------------------------------------------------------------------------------------------------------
# Figure: (a) scatter plot, (b) amplitude divided out and compared with B2(xi)
#-----------------------------------------------------------------------------------------------------------------------
plot_df = result_dataframe[result_dataframe["N"] < 1100]     # first 100 N per strike keeps plot readable
xi_sub = plot_df["xi"].to_numpy()
err_sub = plot_df["scaled_error"].to_numpy()
K_sub = plot_df["K"].to_numpy()

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Panel (a): raw scaled error vs fractional node position
sc0 = ax[0].scatter(xi_sub, err_sub, c=K_sub, cmap="viridis", s=7)
fig.colorbar(sc0, ax=ax[0], label="Strike price K")
ax[0].set_xlabel(r"$\xi$ - fractional node position")
ax[0].set_ylabel(r"$N \times (\mathrm{tree} - \mathrm{BS})$")
ax[0].set_title("(a) Scaled error vs $\\xi$", fontsize=10)

# Panel (b): amplitude divided out, compared to B2(xi)
g_sub = strike_price_payment(S0, K_sub, T, r, sigma)
collapsed = err_sub / (-2 * sigma * np.sqrt(T) * g_sub)

sc1 = ax[1].scatter(xi_sub, collapsed, c=K_sub, cmap="viridis", s=7)
xx = np.linspace(0, 1, 200)
ax[1].plot(xx, second_bernoulli(xx), "k--", lw=1.5, label=r"theory: $B_2(\xi) = \xi^2 - \xi + 1/6$")
ax[1].set_xlabel(r"$\xi$ - fractional node position")
ax[1].set_ylabel(r"$N \times \mathrm{error} / (-2\sigma\sqrt{T} \cdot g(K))$")
ax[1].set_title("(b) Amplitude divided out; offset = smooth $a(K)$", fontsize=10)
ax[1].legend(frameon=False)

fig.tight_layout()
fig.savefig("OE_figure.png", dpi=200)
print("\nSaved: OE_out_of_sample.csv, OE_figure.png")
plt.show()