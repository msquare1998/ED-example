# =========================================================================
#  A simple example for
#       exact diagonalization of the 1D transverse field Ising model
#  Author: Yi-Ming Ding
#  Email: dingyiming@westlake.edu.cn
#  Last updated: May 29, 2024
# =========================================================================*/
import numpy as np
from numpy import kron

pauliX = np.array([[0, 1], [1, 0]])
pauliZ = np.array([[1, 0], [0, -1]])
pauliI = np.eye(2)

def makeLinks_1D(L):
    return np.array([[i, (i + 1) % L] for i in range(L)])

def makeH_1D(L, J, h):
    nQ = L
    H = np.zeros((2 ** nQ, 2 ** nQ))
    links = makeLinks_1D(nQ)
    for idx, ll in enumerate(links):
        hzz = 1
        for q in range(nQ):
            if q in ll:
                hzz = kron(pauliZ, hzz)
            else:
                hzz = kron(pauliI, hzz)
        H += (-J) * hzz

    for i in range(nQ):
        hx = 1
        for q in range(nQ):
            if i == q:
                hx = kron(pauliX, hx)
            else:
                hx = kron(pauliI, hx)
        H += (-h) * hx

    return H

def makeH_1D_shift(L, J, h):
    nQ = L
    H = np.zeros((2 ** nQ, 2 ** nQ))
    links = makeLinks_1D(nQ)
    for idx, ll in enumerate(links):
        hzz = 1
        for q in range(nQ):
            if q in ll:
                hzz = kron(pauliZ, hzz)
            else:
                hzz = kron(pauliI, hzz)
        H += (-J) * (hzz + np.eye(2 ** nQ))

    for i in range(nQ):
        hx = 1
        for q in range(nQ):
            if i == q:
                hx = kron(pauliX, hx)
            else:
                hx = kron(pauliI, hx)
        H += (-h) * (hx + np.eye(2 ** nQ))
    return H

def makeCorr(s0, s1, nQ):
    corr = 1
    for q in range(nQ):
        if s0 != s1:
            if q == s0 or q == s1:
                corr = np.kron(pauliZ, corr)
            else:
                corr = np.kron(np.eye(2), corr)
        else:
            if q == s0:
                corr = np.kron(pauliZ @ pauliZ, corr)
            else:
                corr = np.kron(np.eye(2), corr)
    return corr