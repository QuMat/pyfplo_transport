#!/usr/bin/python
#By Lourdes Amigo and Jorge Facio

import sys
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from update import update_progress
from aux_func import read_data_TDF, read_data_mu, dfdE, integral_n

class transport_tensors():
    """
    Here we have methods that starting from the transport distribution function compute distinct conductivities
    """

    def __init__(self,**kwargs):

        self.root_dir = kwargs.pop('root_dir','./')
        self.prename = kwargs.pop('prename','')

        #build function that interpolates TDF vs T
        self.TDF = {}
        for i in range(3):
            for j in range(3):
                self.ETij,Tij=read_data_TDF(i,j,self.root_dir,self.prename)
                self.TDF[str(i) + str(j)] = interp1d(self.ETij, Tij, kind='linear')


    def compute_sigma(self,T,mu):
        """
        Method that computes the electrical conductivity tensor at a given temperature T and chemical potential mu
        """
        self.sigma = {}
        for i in range(3):
            for j in range(3):

                def integrando(E,T):
                    return -dfdE(E,mu,T)*self.TDF[str(i) + str(j)](E)

                self.sigma[str(i) + str(j)] = quad(integrando, self.ETij[0], self.ETij[-1], args=(T),limit=100)

        return self.sigma
                #tau*electron**2


    def compute_sigma_vs_T_fixed_N(self,Ts = [i * 0.00025 for i in range(1, 140)]):
        """
        Method that computes the electrical conductivity tensor at fixed density as a function of temperature T. 

        Args: 

        Ts: list of temperatures in eV
        
        The data of the chemical potential as a function of temperature is expected in self.root_dir and is built as a functioni in self.mu_int

        Returns:

        It writes the (i,j) component of the conductivity in the file sigma_i_j.dat
        """

        #build function that interpolates mu vs T
        Tt,n,mu = read_data_mu(root_dir=self.root_dir)
        self.mu_int=interp1d(Tt, mu, kind='linear')

        i_T = 0
        for T in Ts:
            mu = self.mu_int(T)
            sigma = self.compute_sigma(T,mu)
            update_progress(i_T * 1. / len(Ts))
            for i in range(3):
                for j in range(3):
                    #print(i,j)
                    name_a = self.root_dir + """/sigma_vs_T_%(i)s_%(j)s.dat"""%locals()
                    if(i_T == 0):
                        print("""#T(K), mu, sigma_%(i)s_%(j)s, error en sigma_%(i)s_%(j)s"""%locals(), file=open(name_a, 'w'))
                    print(T*11604,mu,sigma[str(i)+str(j)][0],sigma[str(i)+str(j)][1],file=open(name_a, 'a'))

    def compute_sigma_vs_mu(self,T,mus = [-0.06 + i * 0.005 for i in range(20)],write_dens=False):
        """
        Method that computes the electrical conductivity tensor as a function of the chemical potential at constant temperature. 

        Args: 

        T: Temperature in eV

        mus: list of chemical potential of interest

        write_dens: boolean, if True, it will add in the output file the density corresponding to each chemical potential.
        

        Returns:

        It writes the (i,j) component of the conductivity in the file sigma_vs_mu_i_j.dat
        """

        i_T = 0
        for mu in mus:
            if(write_dens):
                n = integral_n(mu,T)[0]
                print(n)
            sigma = self.compute_sigma(T,mu)
            update_progress(i_T * 1. / len(mus))
            for i in range(3):
                for j in range(3):
                    name_a = self.root_dir + """/sigma_vs_mu_%(i)s_%(j)s.dat"""%locals()
                    if(i_T == 0):
                        if(write_dens):
                            print("""#T(K), n, mu, sigma_%(i)s_%(j)s, error en sigma_%(i)s_%(j)s"""%locals(), file=open(name_a, 'w'))
                        else:
                            print("""#T(K), mu, sigma_%(i)s_%(j)s, error en sigma_%(i)s_%(j)s"""%locals(), file=open(name_a, 'w'))


                    if(write_dens):
                        print(T*11604,mu,n,sigma[str(i)+str(j)][0],sigma[str(i)+str(j)][1],file=open(name_a, 'a'))
                    else:
                        print(T*11604,mu,sigma[str(i)+str(j)][0],sigma[str(i)+str(j)][1],file=open(name_a, 'a'))

            i_T += 1
			
