import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as spst
from scipy.stats import norm

d1 = lambda S,K,r,sigma,T : (np.log(S/K) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
d2 = lambda S,K,r,sigma,T : (d1(S,K,r,sigma,T) - sigma*np.sqrt(T))


#for newtons method, we need:
callPrice = lambda S,K,r,sigma,T : S*norm.cdf(d1(S,K,r,sigma,T)) - K*np.exp(-r*T)*norm.cdf(d2(S,K,r,sigma,T))
putPrice = lambda S,K,r,sigma,T : K*np.exp(-r*T)*norm.cdf(-d2(S,K,r,sigma,T)) - S*norm.cdf(-d1(S,K,r,sigma,T)) 
vega = lambda S,K,r,sigma,T : S*norm.pdf(d1(S,K,r,sigma,T))*np.sqrt(T)

# newton method:
def newton(optionPrice,S,K,r,T,optionType,initial):
    # start with initial guess
    guess = initial
    
    # Assume the option is a call
    optionFuction = callPrice
    
    # unless it isnt, then dont.
    if optionType == 1:
       optionFuction = putPrice

    # 10 iterations should be enough for mr. newton
    for i in np.arange(10):

        guessChange = (optionPrice-optionFuction(S,K,r,guess,T))/vega(S,K,r,guess,T)
        guess = guessChange + guess
    
    return(guess)

# testing the method
optionprice, S,Kcall,Kput,r,T = 4, 100,110,90,.02,1
sigma =newton(optionprice,S,Kcall,r,T,0,.7)

# this should be 4.
print(callPrice(S,Kcall,r,sigma,T))

# test for puts.
sigma = newton(optionprice,S,Kput,r,T,1,.7)
#should be 4
print(putPrice(S,Kput,r,sigma,T))
