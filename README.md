# README.md
An example of exact diagonalization method for the 1D transverse field Ising model (TFIM) using Python, including calculating

- Energy
- ```lnZ```, where ```Z``` is the partition function
- Spin-spin correlation

This program presents a very brutal way for exact diagonalization, without considering any symmetry, therefore it is simple and can be used for benchmarking other methods like quantum Monte Carlo or density matrix renormalization group.

Two functions are provided
- ```makeH_1D(L, J, h)```: for the standard TFIM
- ```makeH_1D_shift(L, J, h)```: for the shifted TFIM considered in stochastic series expansion

Contact me: dingyiming@westlake.edu.cn
