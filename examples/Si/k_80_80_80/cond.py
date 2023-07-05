import sys,os

from chemical_potential import find_mu
from cond_tensor import transport_tensors
from update import update_progress

def work():

    dirf = '''./'''%locals()
    os.chdir(dirf)

    Ts = [i * 0.001 for i in range(1, 30)]

    if( not os.path.exists('muvsT.dat')):
        i_T = 0
        for T in Ts:
            update_progress(i_T * 1./len(Ts))
            mu = find_mu(T,Nobj=8,mu0=-1,mu1=1)
            i_T += 1
    my_transport = transport_tensors()
    my_transport.compute_sigma_vs_T_fixed_N(Ts=Ts)
    my_transport.compute_sigma_vs_mu(T=300./11604,mus=[-2 + i*0.02 for i in range(200)])

if __name__ == '__main__':

    work()

