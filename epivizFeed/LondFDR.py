'''
based on https://github.com/dsrobertson/onlineFDR
'''

import numpy as np
import scipy
import math

class LondFDR(object):
    def __init__(self, R=0, N=1):
        self.R = R
        self.N = N

    def batchLondStar(self, pval, alpha = 0.05, betai = None):
        # print("######## IN LOND ##########")
        # print(pval)
        
        N = len(pval)
        alphai = [None for x in range(0,N)]
        R = [0 for x in range(0,N)]

        D = max(self.R, 1)
        for i in range(0, N): 
            alphai[i] = (0.07720838*alpha*math.log(max(self.N,2)) / (self.N*math.exp(math.sqrt(math.log(self.N))))) * D
            self.N += 1
            R[i] = float(pval[i]) <= alphai[i]

        self.R += np.sum(R)
        return (pval, alphai, R)

    def lond(self, pval, alpha = 0.05, betai = None, dep=False, pOrderRand = True, algoOrig = True):
        
        N = len(pval)
        seq_N = np.arange(1, N+1, 1)

        if betai is None:
            betai = []
            for x in range(1, N+1):
                tbetai = 0.07720838*alpha*math.log(max(x,2)) / (x*math.exp(math.sqrt(math.log(x))))
                betai.append(tbetai)
        elif min(betai) < 0:
            raise Exception("BetaiCannotBeNegative")
        elif sum(betai) > alpha:
            raise Exception("SumGreaterThanAlpha")

        den = None

        if dep:
            den = np.cumsum(1/seq_N)
            betai = betai / den

        ### LOND Procedure

        R = [0 for x in range(0,N)]
        alphai = [0 for x in range(0,N)]

        alphai[0] = betai[0]
        R[0] = pval[0] <= alphai[0]

        if N == 1:
            return (pval, alphai, R)

        D = R[0]

        if not algoOrig:
            for x in range(1, N):
                alphai[x] = betai[x]*max(D,1)
                R[x] = pval[x] <= alphai[x]
                D = D + R[x]
        else:
            for x in range(1, N):
                alphai[x] = betai[x]*(D+1)
                R[x] = pval[x] <= alphai[x]
                D = D + R[x]

        return (pval, alphai, R)

    def londStar(self, pval, alpha = 0.05, betai = None, version = "dep", decision = None, lags = None, batchSize = None):

        N = len(pval)
        seq_N = np.arange(1, N+1, 1)

        if betai is None:
            betai = []
            for x in range(1, N+1):
                tbetai = 0.07720838*alpha*math.log(max(x,2)) / (x)*math.exp(math.sqrt(math.log(x)))
                betai.append(tbetai)
        elif min(betai) < 0:
            raise Exception("BetaiCannotBeNegative")
        elif sum(betai) > alpha:
            raise Exception("SumGreaterThanAlpha")

        L = lags
        
        R = alphai = np.repeat(0, N)
        
        alphai[0] = betai[0]
        R[0] = pval[0] <= alphai[0]
        
        if N == 1:
            return (pval, lags, alphai, R)
        
        for x in range(1, N):
            D = max(sum(R[np.arange(1, x, 1)] == 1 & np.arange(1, x, 1) <= x-1-L[x]), 1)
            
            alphai[x] = betai[x]*D
            R[x] = pval[x] <= alphai[x]
        
        return (pval, lags, alphai, R)