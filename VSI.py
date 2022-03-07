import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as spst
from scipy.stats import norm
from scipy.optimize import minimize 
import yfinance as yf

#get data from yfinance
aapl= yf.Ticker("aapl")
opt = aapl.option_chain('2022-08-19')

#only grab out of the money options
calls = opt.calls.loc[opt.calls.inTheMoney == False]
puts = opt.puts.loc[opt.puts.inTheMoney == False]

#plot
strikesCalls = calls.strike
impliedVolsCalls = calls.impliedVolatility
strikesPuts = puts.strike
impliedVolsPuts = puts.impliedVolatility

strikes = pd.concat([strikesCalls,strikesPuts],axis=0)
impliedVols = pd.concat([impliedVolsCalls,impliedVolsPuts],axis=0)

#define our functions to mimimize
sviUnvectorized = lambda k,a,b,rho,m,sigma : a+b*(rho*(k-m)+np.sqrt((k-m)**2 + sigma**2))

svi = np.vectorize(sviUnvectorized)

def funcUnvectorized(parameters):
    var = impliedVols
    k = strikes
    cost = sum((var - svi(k,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]))**2)
    return(cost)
# instead of using a minimize routine for the initial guesses,
# initial guesses are calculated based on section 2 of the paper "Initial Guesses for SVI Calibration" by Le Floc'h
# Essentially, the asymptotes are approximated and parameters are based on that.

# These are my asymptotic approximations: ar and al are the slopes of the right and left asymptotes; br and bl are the intercepts.
#I calculated these by hand; could be automated by calculating a line of best fit for the 5 rightmost and then for the 5 leftmost data points.
# right asymptote = ar(k) + br = -(.069/20)*k+.83
# left asymptote = al(k) + bl =(.021/30)*kvals+.117


# these are the parameters in terms of the asymptotic approximations (Le Floc'h):
# b = .5(ar - al)
# rho = (ar+al)/(ar-al)
# a = bl - al((bl-br)/(al-ar))
# m = (-bl + a)/(b(1-rho))

# sigma needs the minimum implied variance, minVar
# sigma = (minvar - a)/(b(sqrt(1-rho^2)))

#manually calculated the guesses below; obviously this could be automated.
minimum = minimize(funcUnvectorized,x0=[.237265,-.002075,.668,-506,3.86],bounds=[(None,None),(.01,None),(-.999999,1),(None,None),(0.005,None)])
print(minimum)

#-.002075

#sometimes the minimizing function likes returning a straight line. Here's some parameter values from a successful return, frozen here.
parametersFrozen = np.array([ -12.40526196,    2.6829,   -0.9966945, -505.440226, 58.38874569])

kvals=np.linspace(75,265,1500)
plt.plot(strikes,impliedVols,"x", label = 'actual')
plt.plot(kvals,svi(kvals,parametersFrozen[0],parametersFrozen[1],parametersFrozen[2],parametersFrozen[3],parametersFrozen[4]), label='Paramters that Work')
plt.plot(kvals,svi(kvals,minimum.x[0],minimum.x[1],minimum.x[2],minimum.x[3],minimum.x[4]), label = 'Unreliable Automatic SVI Parameterization')
plt.title("options that expire on 2022-08-19")
plt.xlabel("strikes")
plt.ylabel("implied vol")
plt.legend()
plt.show()

#a,b,rho,m,sigma
