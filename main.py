# =========================================================================
#  A simple example for
#       exact diagonalization of the 1D transverse field Ising model
#  Author: Yi-Ming Ding
#  Email: dingyiming@westlake.edu.cn
#  Last updated: May 29, 2024
# =========================================================================*/
import numpy as np
from numpy import trace as tr
from numpy import exp, log
from scipy.linalg import expm
import aux
from scipy import sparse

nQ = L = 6
beta = 5
J = 1
h = 1

if __name__ == "__main__":
    print(f"# L = {L}, beta = {beta}, J = {J}, h = {h}")
    H = aux.makeH_1D(L, J, h)
    """ -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        Calculating the energy
    -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ """
    print(f"# E = ", end="")
    energy = tr(expm(-beta * H) @ H) / tr(expm(-beta * H))
    assert energy.imag < 1e-9, f"Imaginary number occurs when calculating energy, with energy.imag = {energy.imag}"
    print(f"{energy.real}")

    """ -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        Calculating lnZ
    -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ """
    print(f"# lnZ = ", end="")
    H_sp = sparse.csr_matrix(H)
    val_min, _ = sparse.linalg.eigsh(H_sp, k=2 ** (nQ - 1), which="SA")
    val_max, _ = sparse.linalg.eigsh(H_sp, k=2 ** (nQ - 1), which="LA")
    val = np.concatenate((val_min, np.sort(val_max)), axis=0)
    Z = np.array([exp(-beta * x) for x in val]).sum()
    print(f"{log(Z)}")

    """ -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        Calculating the spin-spin correlation
    -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ """
    print("# Correlations:", end="")
    for i in range(L):
        corr_op = aux.makeCorr(0, i, L)
        corr = tr(corr_op @ expm(-beta * H)) / tr(expm(-beta * H))
        if corr.imag < 1e-10:
            corr = corr.real
        else:
            raise Exception("\nImaginary errors.")
        if i % (int(L / 2)) == 0:
            print(f"\n\tcorr(0, %d) = %.7f" % (i, corr), end="\t\t\t")
        else:
            print(f"\tcorr(0, %d) = %.7f" % (i, corr), end="\t\t\t")
    print("")