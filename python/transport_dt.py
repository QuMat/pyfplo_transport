#!/usr/bin/python
#By Lourdes Amigo and Jorge facio

import sys,os
import numpy as np
import numpy.linalg as LA
import h5py
import pyfplo.slabify as sla
from update import update_progress

import time

vel_to_angs_eV = 0.529177249
h_eVfs = 4.135667696 #eV fs
from_G0 = 7.748091729 #10e-5 S

class bz_conductivity():

    # Here we have methods that implement the calculation of the transport distribution tensor, similarly than done in Computer Physics Communications 185 (2014)
    def __init__(self,**kwargs):
        """
        Method that initializes our transport calculation. 

        We assume that the pyfplo calculation was done including the spin orbit coupling.

        Args:

        slabify: pyfplo object built from a Wannier Hamiltonian
        """
        self.S = kwargs['slabify']
        self.ms  = 0 #we always have spin-orbit coupling
        self.gauge = kwargs.pop('gauge','periodic')

        self.origin = kwargs.pop('origin',[0.,0.,0.])
        self.root = kwargs.pop('root_name','.')
        R=self.S.hamdataCell()
        G=self.S.hamdataRCell()
        self.G1=G[:,0]
        self.G2=G[:,1]
        self.G3=G[:,2]
        self.dG1=LA.norm(self.G1)
        self.dG2=LA.norm(self.G2)
        self.dG3=LA.norm(self.G3)
        
        self.real_VOL_fplo = np.abs(np.linalg.det(R))
        self.real_VOL_Ang3 = self.real_VOL_fplo * 0.529177249**3 

        self.write_info()
        #print("real space unit cell: ", R)
        #print("real space unit cell volume in fplo^3: ", self.real_VOL_fplo)
        print("real space unit cell volume in Angs^3: ", self.real_VOL_Ang3)
        #print "scale: ", self.S.kscale


    def compute_E_and_velocities(self,k,test_diag=False):
        """
        Computes eigenvalues and velocities
        
        Args:

        k: np array with the k-points coordinate in the units used in the +hamdata file

        Returns:

        Output of sla.diagonalize(makef=True) so it is in units of eV Bohr 
        """
        (Hk,dHk) = self.S.hamAtKPoint(k,self.ms,gauge=self.gauge,makedhk=True)
        (E,U) = self.S.diagonalize(Hk)

        Udag = np.conjugate(np.transpose(U))

        if(test_diag):
            Ed = LA.multi_dot([Udag, Hk, U])
            assert(LA.norm(Ed-np.diag(np.diagonal(Ed)))<1e10)

        vx = np.diagonal(LA.multi_dot([Udag, dHk[0], U]))
        vy = np.diagonal(LA.multi_dot([Udag, dHk[1], U]))
        vz = np.diagonal(LA.multi_dot([Udag, dHk[2], U]))

        return E,[vx,vy,vz]

    def write_info(self):
        """
        Method to write some global information in our ar.hdf5 file.
        """
        root = self.root
        self.ar = h5py.File("""%(root)s/ar.hdf5"""%locals(),"w")
        run_info = self.ar.create_group("run_info")
        run_info.create_dataset('kscale',data=self.S.kscale)
        self.ar.close()

    def build_box(self,k_mesh):
        """
        Build a 3d mesh of the BZ
        """
        box=sla.BoxMesh()
        box.setBox(xaxis=self.G1,yaxis=self.G2,zaxis=self.G3,origin=self.origin)

        box.setMesh(nx=k_mesh[0],xinterval=[-0.5*self.dG1,(0.5-1./k_mesh[0])*self.dG1], ny=k_mesh[1],yinterval=[-0.5*self.dG2,(0.5-1./k_mesh[1])*self.dG2], nz=k_mesh[2],zinterval=[-0.5*self.dG2,(0.5-1./k_mesh[2])*self.dG3])

        kpoints = box.mesh(self.S.kscale) #after this the points are in Bohr
        return kpoints


    def tdf(self,k_mesh=[10,10,1],e_mesh=[-0.1 + 0.01*i for i in range(20)],out_folder='.'):
        """
        Method to compute the transport distribution function defined in Eq. 8 of Ref. 1.
        The relaxation time is set to 1 fs, so it should best in this units in following postprocessing based on the obtained TDF.

        Args:
        
        k_mesh: list of the three integers describding the discretization of the BZ.
        e_mesh: energy mesh at which the TDF will be computed

        Returns:

        A list of 3x3 np.matrix objects containing for each energy in e_mesh the TDF tensor.
        """
      
        Sigma = [np.matrix(np.zeros((3, 3)),dtype=complex) for e in e_mesh] #one 3x3 matrix for each energy in e_mesh
        DOS = [ 0 for e in e_mesh] #one 3x3 matrix for each energy in e_mesh

        dE = e_mesh[1] - e_mesh[0]

        factor_dos = 1./(k_mesh[0]*k_mesh[1]*k_mesh[2]) / dE
        factor_sigma = (2. * np.pi**2 * 1e-3 * from_G0 / h_eVfs) * factor_dos  * vel_to_angs_eV**2   / self.real_VOL_Ang3 

        kpoints = self.build_box(k_mesh) #in Bohr -1
        i_k = 0
        for k in kpoints:
            update_progress(i_k*1./len(kpoints))
            E,V = self.compute_E_and_velocities(k)

            for n in range(len(E)):
                if(E[n] >= e_mesh[0] and E[n] <= e_mesh[-1]):

                    energy_index = np.int(np.floor((E[n] - e_mesh[0])/(e_mesh[1]-e_mesh[0])))
                    DOS[energy_index] += 1 * factor_dos #with this factor the DOS should integrate to the number of electrons.

                    for i in range(3):
                        for j in range(3):
                            Sigma[energy_index][i,j] += V[i][n] * V[j][n] * factor_sigma

            i_k += 1
        for i in range(3):
            for j in range(3):
                file = open("""%(out_folder)s/TDF_%(i)s_%(j)s.dat"""%locals(),'w')
                for n in range(len(e_mesh)):
                    file.write(str(e_mesh[n]) + ' ' + str(Sigma[n][i,j].real) + ' ' + str(Sigma[n][i,j].imag) + '\n')
                file.close()


        occupancy = 0.0
        file = open("""%(out_folder)s/dos.dat"""%locals(),'w')
        for n in range(len(e_mesh)):
            if(e_mesh[n] <= 0):
                occupancy += DOS[n] * dE
            file.write(str(e_mesh[n]) + ' ' + str(DOS[n]) + '\n')

        print("The ocupancy|density are: ", occupancy , '|' ,occupancy)

        file.close()
        root = self.root
        self.ar = h5py.File("""%(out_folder)s/%(root)s/ar.hdf5"""%locals(),"a")
        G_data = self.ar.create_group("transport")
        G_data.create_dataset('tdf',data=Sigma)

        return Sigma

    def tdf_ioffe_regel(self,lat_params,k_mesh=[10,10,1],e_mesh=[-0.1 + 0.01*i for i in range(20)],out_folder='.'):
        """
        Computes the Ioffe Regel transport distribution function, wich is defined by assuming the relaxation time to be tau = lat_param/v_f 

        Args:

            lat_params: It should be a list with three numbers associated with the lattice parameters in Angstroms.
        """

        Sigma = [np.matrix(np.zeros((3, 3)),dtype=complex) for e in e_mesh] #one 3x3 matrix for each energy in e_mesh
        DOS = [ 0 for e in e_mesh] #one 3x3 matrix for each energy in e_mesh

        print('lat_params: ', lat_params)

        dE = e_mesh[1] - e_mesh[0]
        factor_dos = 1./(k_mesh[0]*k_mesh[1]*k_mesh[2]) / dE
        factor_sigma = (2. * np.pi**2 * 1e-3 * from_G0 / h_eVfs) * factor_dos  * vel_to_angs_eV**2   / self.real_VOL_Ang3 

        kpoints = self.build_box(k_mesh) #in Bohr -1
        i_k = 0
        for k in kpoints:
            update_progress(i_k*1./len(kpoints))
            E,V = self.compute_E_and_velocities(k)

            for n in range(len(E)):
                if(E[n] >= e_mesh[0] and E[n] <= e_mesh[-1]):

                    energy_index = np.int(np.floor((E[n] - e_mesh[0])/(e_mesh[1]-e_mesh[0])))
                    DOS[energy_index] += 1 * factor_dos #with this factor the DOS should integrate to the number of electrons.

                    for i in range(3):
                        for j in range(3):
                            average_fermi_velocity = np.sqrt(V[i][n]**2 + V[j][n]**2 )/ np.sqrt(2.) * vel_to_angs_eV / (h_eVfs / (2* np.pi))
                            
                            average_lat_param = 0.5 * (lat_params[i] + lat_params[j])

                            #here the idea is that if the velocity is very small, we can put tau to e.g. 1 because this state will not contribute to the conductivity and we would like not have a diverging tau involved.
                            if(average_fermi_velocity > 1e-8):
                                tau = average_lat_param / average_fermi_velocity
                            else:
                                tau = 1

                            Sigma[energy_index][i,j] += tau * V[i][n] * V[j][n] * factor_sigma

            i_k += 1
        for i in range(3):
            for j in range(3):
                file = open("""%(out_folder)s/IoffeRegel_TDF_%(i)s_%(j)s.dat"""%locals(),'w')
                for n in range(len(e_mesh)):
                    file.write(str(e_mesh[n]) + ' ' + str(Sigma[n][i,j].real) + ' ' + str(Sigma[n][i,j].imag) + '\n')
                file.close()


        occupancy = 0.0
        file = open("""%(out_folder)s/dos.dat"""%locals(),'w')
        for n in range(len(e_mesh)):
            if(e_mesh[n] <= 0):
                occupancy += DOS[n] * dE
            file.write(str(e_mesh[n]) + ' ' + str(DOS[n]) + '\n')

        print("The ocupancy|density are: ", occupancy , '|' ,occupancy)

        file.close()
        root = self.root
        self.ar = h5py.File("""%(out_folder)s/%(root)s/ar.hdf5"""%locals(),"a")
        G_data = self.ar.create_group("transport")
        G_data.create_dataset('tdf',data=Sigma)

        return Sigma
