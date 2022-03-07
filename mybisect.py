import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as spst
from scipy.stats import norm


d1 = lambda S,K,r,sigma,T : (np.log(S/K) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
d2 = lambda S,K,r,sigma,T : (d1(S,K,r,sigma,T) - sigma*np.sqrt(T))


#for newtons method, we need:
callPrice = lambda S,K,r,sigma,T : S*norm.cdf(d1(S,K,r,sigma,T)) - K*np.exp(-r*T)*norm.cdf(d2(S,K,r,sigma,T))
putPrice = lambda S,K,r,sigma,T : K*np.exp(-r*T)*norm.cdf(-d2(S,K,r,sigma,T)) - S*norm.cdf(-d1(S,K,r,sigma,T)) 

def bisection(optionPrice,S,K,r,T,optionType,upperBoundGuess):
    a=0
    b=upperBoundGuess
    optionFunction = callPrice
    if optionType == 1 :
        optionFunction = putPrice

    for i in np.arange(100):
        midpoint = (a+b)/2
        if optionPrice<optionFunction(S,K,r,midpoint,T):
            b = midpoint
        elif optionPrice>optionFunction(S,K,r,midpoint,T):
            a = midpoint
        elif optionPrice == optionFunction(S,K,r,midpoint,T):
            print('done early')
            return(midpoint)
        else:
            print('error')
    return(midpoint)


# Testing: Call
sigma = bisection(4,100,110,.02,1,0,10)
print("sigma is " ,sigma)

# this should be 4
print(callPrice(100,110,.02,sigma,1))

# Testing: Put
sigma = bisection(4,100,90,.02,1,1,10)
print("sigma is" , sigma)

# this should be 4
print(putPrice(100,90,.02,sigma,1))

