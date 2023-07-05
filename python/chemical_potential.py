#!/usr/bin/python
#By Lourdes Amigo & Jorge Facio.

import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d

from aux_func import interpola_n
	
def find_mu(T,Nobj=12,root_dir='.',mu0 = -0.4, mu1 = 0.37,Nmu=100,method='bisect'):
    """
    Method that finds chemical potential starting from data of the density as a function of the chemical potential

    Args:

    T: temperature in eV
    Nobj: target density 
    root_dir: directory where to find the file nvsmu.dat which contains the density as a function of the chemical potential
    mu0, mu1: starting points fo the chemical potential for the bisection method to find roots.
    Nmu: number of points for the discretization of the chemical potential domain.

    Returns:

    The chemical potential found, which is also appended to to a file called muvsT.dat
    """
    name_a = """%(root_dir)s/nvsmu.dat"""%locals()
    name_c = """%(root_dir)s/muvsT.dat"""%locals()
    delta_mu = (mu1-mu0)/Nmu
    mus = [mu0*1.001 + i * delta_mu for i in range(Nmu+2)]
    ny = []
    for mu in mus:
        nint=integral_n(mu,T,root_dir)[0]
        print(T*11604,mu,nint, file=open(name_a, 'a'))
        ny.append(nint)
    
    nvsmu=interp1d(mus, ny, kind='linear')


    def interpola_n(mu,nn):
        return nvsmu(mu)-nn

    solucion = root_scalar(interpola_n,args=(Nobj), bracket=[mu0,mu1], method='bisect')
    print(T*11604,Nobj,solucion.root,file=open(name_c,'a'))

    return solucion.root
