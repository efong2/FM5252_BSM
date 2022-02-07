import numpy as np
from scipy.stats import norm

# Functions necessary include call price, put price, call delta, put delta,
# gamma, vega, call theta, put theta, call rho, put rho

# d1 and d2 as handy dandy lambda functions
d1 = lambda S,K,r,sigma,T : (np.log(S/K) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
d2 = lambda S,K,r,sigma,T : (d1(S,K,r,sigma,T) - sigma*np.sqrt(T))

# functions below take a stock price, S; a strike, K; a rate, r; a volatility, sigma; and a maturity T.

# returns callprice according to BSM. 
callPrice = lambda S,K,r,sigma,T : S*norm.cdf(d1(S,K,r,sigma,T)) - K*np.exp(-r*T)*norm.cdf(d2(S,K,r,sigma,T))

# returns put price according to BSM
putPrice = lambda S,K,r,sigma,T : K*np.exp(-r*T)*norm.cdf(-d2(S,K,r,sigma,T)) - S*norm.cdf(-d1(S,K,r,sigma,T)) 

# returns call delta, the change in option price given a change in the underlying price
callDelta = lambda S,K,r,sigma,T : norm.cdf(d1(S,K,r,sigma,T))

# returns put delta, the change in option price given a change in the underlying price
putDelta = lambda S,K,r,sigma,T : norm.cdf(d1(S,K,r,sigma,T))-1

# returns gamma, the change in delta given a change in the underlying
gamma = lambda S,K,r,sigma,T : norm.pdf(d1(S,K,r,sigma,T))/(S*sigma*np.sqrt(T))

# returns vega, the change in the option price given a change in volatility
vega = lambda S,K,r,sigma,T : S*norm.pdf(d1(S,K,r,sigma,T))*np.sqrt(T)

# returns Theta, the change in price of a put or call option as time passes. different for calls and puts
callTheta = lambda S,K,r,sigma,T : -S*norm.pdf(d1(S,K,r,sigma,T))*sigma/(2*np.sqrt(T)) - r*K*np.exp(-r*T)*norm.cdf(d2(S,K,r,sigma,T))
putTheta = lambda S,K,r,sigma,T:-S*norm.pdf(d1(S,K,r,sigma,T))*sigma/(2*np.sqrt(T)) + r*K*np.exp(-r*T)*norm.cdf(-d2(S,K,r,sigma,T))

# returns rho, the change in price of an option with respect to a change in interest rates.
callRho = lambda S,K,r,sigma,T : K*T*np.exp(-r*T)*norm.cdf(d2(S,K,r,sigma,T))
putRho = lambda S,K,r,sigma,T : -K*T*np.exp(-r*T)*norm.cdf(-d2(S,K,r,sigma,T))