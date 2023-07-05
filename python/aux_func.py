#!/usr/bin/python
#By Lourdes Amigo and Jorge Facio

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def fermi(E,mu,T):
    """
    Equilibrium Fermi distribution.

    Args:

    E: energy
    mu: chemical potential
    T: Temperature in eV 
    """
    return 1. / (1+np.exp(np.float128((E-mu)/T)))
	
def dfdE(E,mu,T):
    """
    Derivative of Fermi distribution function, with everything in eV.
    """
    return -1*np.exp((E-mu)/T)/(T*(np.exp((E-mu)/T)+1)**2)

def integral_n(mu,T,root_dir="./"):
    """
    Method that integrates the density of states times the equilibrium Fermi function.

    Args:

    mu: chemical potential in eV
    T: Temperature in eV
    root_dir: directory where to find the density of states.
    
    """
    E,y = read_dos(root_dir=root_dir)
    DOS = interp1d(E, y, kind='linear')
	
    def integrando(E,mu,T):
        return fermi(E,mu,T)*DOS(E)

    return quad(integrando, E[0], E[-1], args=(mu, T),limit=100)

def read_dos(root_dir="./"):
    """
    Method that reads the density of states.

    Args:

    root_dir: directory where to find the file dos.dat

    Return:

    x,y: list of energies and DOS.
    """
    file = open("""%(root_dir)s/dos.dat"""%locals(),'r')
    data = file.readlines()
    file.close()

    x  = []
    y = []
    for i in range(0,len(data),1):
        line = data[i]
        a = line.split()
        x.append(eval(a[0]))
        y.append(eval(a[1]))

    return x,y
	



def read_data_TDF(i,j,root_dir=".",prename=''):
    """
    Method that reads a the (i,j) component of the TDF file which should be present in the root_dir directory
    """
    file = open("""%(root_dir)s/%(prename)sTDF_%(i)s_%(j)s.dat"""%locals(),'r')
    data2 = file.readlines()
    file.close()

    ETij  = []
    Tij = []
    for i in range(1,len(data2)):
        vals = list(map(eval,data2[i].split()))
        ETij.append(vals[0])
        Tij.append(vals[1])

    return ETij,Tij
	
	
def read_data_mu(root_dir="."):
    """
    Method that reads the chemical potential as a function of temperature, which should be given in the file muvsT.dat in the root_dir directory
    """
    file = open("""%(root_dir)s/muvsT.dat"""%locals(),'r')
    data3 = file.readlines()
    file.close()

    T  = []
    n = []
    mu = []
    for i in range(0,len(data3)):
        vals = list(map(eval,data3[i].split()))
        print(vals)
        fix_units = 1
        if(vals[0] >1):
            fix_units = 11604
        T.append(vals[0]/fix_units)
        n.append(vals[1])
        mu.append(vals[2])
    return T,n,mu
	
def dfdE(E,mu,T):
    """
    Derivative of Fermi distribution function, with everything in eV.
    """
    return -1*np.exp((E-mu)/T)/(T*(np.exp((E-mu)/T)+1)**2)



def read_data_TDF(i,j,root_dir=".",prename=''):
    """
    Method that reads a the (i,j) component of the TDF file which should be present in the root_dir directory
    """
    file = open("""%(root_dir)s/%(prename)sTDF_%(i)s_%(j)s.dat"""%locals(),'r')
    data2 = file.readlines()
    file.close()

    ETij  = []
    Tij = []
    for i in range(1,len(data2)):
        vals = list(map(eval,data2[i].split()))
        ETij.append(vals[0])
        Tij.append(vals[1])

    return ETij,Tij
	
	
def read_data_mu(root_dir="."):
    """
    Method that reads the chemical potential as a function of temperature, which should be given in the file muvsT.dat in the root_dir directory
    """
    file = open("""%(root_dir)s/muvsT.dat"""%locals(),'r')
    data3 = file.readlines()
    file.close()

    T  = []
    n = []
    mu = []
    for i in range(0,len(data3)):
        vals = list(map(eval,data3[i].split()))
        print(vals)
        fix_units = 1
        if(vals[0] >1):
            fix_units = 11604
        T.append(vals[0]/fix_units)
        n.append(vals[1])
        mu.append(vals[2])
    return T,n,mu
	



