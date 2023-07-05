import sys,os


import pyfplo.slabify as sla
from transport_dt import bz_conductivity



def work():

    hamdata='+hamdata'
    s=sla.Slabify()
    s.object='3d'
    s.printStructureSettings()
    s.prepare(hamdata)

    Nx =80  
    Ny =80
    Nz =80

    k_mesh=[Nx,Ny,Nz]
    e_step = 0.05
    e_0 = -13
    e_f = 5
    n_e = int((e_f-e_0)/e_step)
    e_mesh=[e_0 + e_step*i for i in range(n_e)]
    dirf = '''k_%(Nx)s_%(Ny)s_%(Nz)s'''%locals()
    os.mkdir(dirf)

    fe = bz_conductivity(slabify=s)
    fe.tdf(k_mesh = k_mesh, e_mesh = e_mesh, out_folder = dirf)



if __name__ == '__main__':

    work()

