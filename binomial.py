import numpy as np
import scipy.stats as spst

# note, optionsMatrix must have calculated values for its last column!

# function takes in an all but empty matrix, starting coordinates, a strike K, a rate r, a timestep t, a 0 or 1 for nationality,
# a 0 or 1 for call_or_put, and a matrix of underlying prices.

# function returns a solved matrix of options prices.
def optionEvaluator(optionsMatrix,coordinates,K,r,t,p,nationality,call_or_put,underlying):
    coordinatesUpper = (coordinates[0],coordinates[1]+1)
    coordinatesLower = (coordinates[0]+1,coordinates[1]+1)
    # check to see if we're on the tree. if we're off-tree, call the function at the top of the previous row.
    if coordinates[0]>coordinates[1]:
        coordinates = (0,coordinates[1])
        optionEvaluator(optionsMatrix,coordinates,K,r,t,p,nationality,call_or_put,underlying)
    # check to see if we're standing on a nan.
    elif np.isnan(optionsMatrix[coordinates]):
        # if we are, check to see if lower coord is nan
        if np.isnan(optionsMatrix[coordinatesLower]):
            # if the lower coord is a nan, move to the lower coord and call the function again.
            optionEvaluator(optionsMatrix,coordinatesLower,K,r,t,p,nationality,call_or_put,underlying)
        # if lower coord is not a nan, calculate. then, step back and call the function again.
        else:
            expected = p*(optionsMatrix[coordinatesUpper]) + (1-p)*(optionsMatrix[coordinatesLower])
            discount = np.exp(-r*t)
            optionsMatrix[coordinates] = expected*discount
            
            # if evaluating an american option, compare the exercised value against the calculated value.
            if nationality==1:
                optionsMatrix[coordinates] = np.maximum(optionsMatrix[coordinates],(underlying[coordinates]-K)*(-1)**call_or_put)
            
            # if 0,0 is still nan, step back and call the function.
            if np.isnan(optionsMatrix[0,0]):
                coordinates = (coordinates[0],coordinates[1]-1)
                optionEvaluator(optionsMatrix,coordinates,K,r,t,p,nationality,call_or_put,underlying)
        return(optionsMatrix)

    

# nationality: 0=euro, 1=american
# call_or_put: 0=call, 1=put

def binomialModel(S,K,r,sigma,T,N,nationality,call_or_put) : 
    t = T/N
    u = np.exp(sigma*np.sqrt(t))
    d = 1/u
    p = (np.exp(r*t)-d)/(u-d)

    # underlying tree
    uPowers = np.zeros((N+1,N+1))
    uPowers = uPowers + np.arange(N+1)
    uPowers = (uPowers.transpose() - np.arange(N+1)*2).transpose()
    uPowers = np.triu(uPowers)
    underlying = S*u**uPowers

    #Option Tree
    optionsMatrix = np.empty((N+1,N+1))
    optionsMatrix[:] = np.nan
    # Calculate final column of option prices
    optionsMatrix[:,-1] = np.maximum(0,(-1)**call_or_put*(underlying[:,-1]-K))


    return(optionEvaluator(optionsMatrix,(0,N-1),K,r,t,p,nationality,call_or_put,underlying))

# European call over 4 time steps.
print('European Call')
print(binomialModel(100,105,.04,.2,1,4,0,0))

# American call over 4 time steps.
print('American Call')
print(binomialModel(100,105,.04,.2,1,4,1,0))

# European put over 4 time steps.
print('European Put')
print(binomialModel(100,105,.04,.2,1,4,0,1))

# American put over 4 time steps
print('American Put')
print(binomialModel(100,105,.04,.2,1,4,0,1))
