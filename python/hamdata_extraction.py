from __future__ import division
import numpy as np
import pickle
import sys
import pylab as py

from update import update_progress

py.ion()


# HAMDATA FILE PROCESSING
# This is a modification of a code written by I.C. Fulga

known_data_types = ['from_fplo',]


def process_hamdata(infile='.', outfile='.', ftype='from_fplo',
                    return_data=False, keep_only_atom_name=True):
    """ the ftype string is supposed to make this function work with a variety
    of different ab-initio codes, which maybe encode the Wannier Hamiltonian
    in different ways...
    
    Returned file format is: a, atom_pos, atom_names, atom_orbs, total_hamdata
        --- 3 lattice vectors, a_1, a_2, a_3, arranged as rows of a 3x3 matrix
        --- list of 3 vectors indicating the atom positions in the unit cell
        --- list of names of the atoms (strings)
        --- list of integers. indicates the number of orbitals per atom
        --- dictionary with keys: 'atom_i atom_j integer_distance_ij'
                    and values: hamiltonian_submatrix.
                This describes a hopping FROM atom_j TO atom_i, which are
                separated by integer_distance_ij unit cells. For example,
                2 unit cells in the a_1 direction, 3 unit cells in the a_2
                direction, and -5 unit cells in the a_3 direction. Zero total
                distance means onsite term. The hamiltonian_submatrix is the
                onsite/hopping term connecting those particular atoms. It is a
                rectangular matrix if the two atoms have different nr. orbitals.
    """
    if ftype not in known_data_types:
        raise ValueError('Unrecognized data type.')

    if ftype == 'from_fplo':
        with open(infile, 'r+') as f:
            lines = f.readlines()
            i = 0
            while not lines[i].startswith('lattice_vectors'):
                i += 1

            # extract 3 lattice vectors
            a = np.zeros((3, 3), dtype=float)
            for j in range(3):
                a[j,:] = [float(s) for s in lines[i+j+1].split()]

            inva = np.linalg.inv(a.T) # used to determine 'integer_distance_ij'

            while not lines[i].startswith('nwan'):
                i += 1

            nwan = int(lines[i+1]) # total number of Wannier orbitals

            while not lines[i].startswith('wannames'):
                i += 1

            wannames = [] # names of Wannier orbitals
            for ind in range(nwan):
                if lines[i+1+ind][-1:] == '\n':
                    wannames.append(lines[i+1+ind][:-1])
                else:
                    wannames.append(lines[i+1+ind])

            while not lines[i].startswith('wancenters'):
                i += 1

            wcen = np.loadtxt(lines[i+1:i+1+nwan]) # positions of all orbitals

            # remove the spin component and orbital name from Wannier names
            if keep_only_atom_name:
                for ind in range(nwan):
                    wannames[ind] = wannames[ind][:wannames[ind].find(' ')]

            atom_pos = [] # atom positions
            atom_names = [] # atom names
            atom_orbs = [] # nr. orbitals for each atom
            orb_ind = [] # to which atom does each orbital belong
            orb_ind_start = [] # starting orbital index of each atom
            orb_ind_stop = [] # final orbital index of each atom

            orb_ind_start.append(0) # first atom starts with orbital nr. 0

            atom_pos.append(wcen[0, :])
            atom_names.append(wannames[0])
            orb_current = wcen[0, :]
            orb_counter = 0
            nratoms = 1
            for ind in range(nwan):
                orb_counter += 1
                if np.linalg.norm(np.abs(wcen[ind, :] - orb_current)) >= 1e-5:
                    atom_pos.append(wcen[ind, :])
                    atom_names.append(wannames[ind])
                    orb_current = wcen[ind, :]
                    nratoms +=1
                    atom_orbs.append(orb_counter-1)
                    orb_counter = 1
                    orb_ind_start.append(ind)
                    orb_ind_stop.append(ind)

                orb_ind.append(nratoms - 1)

            # last atom
            atom_orbs.append(orb_counter)
            orb_ind_stop.append(nwan)

            # determine which parts of the file contain the Wannier data
            data_start = []
            data_stop = []
            for ind in range(len(lines)):
                if lines[ind].startswith('Tij, Hij:'):
                    data_start.append(ind)

                if lines[ind].startswith('end Tij, Hij:'):
                    data_stop.append(ind)

            if len(data_start) != nwan**2:
                raise ValueError('The number of entries in the data file' +
                                 ' does not mach the number of' +
                                 ' Wannier orbitals.')

            # begin extracting Hamiltonain elements (hamiltonian_submatrix)
            total_hamdata = {}

            datalen = len(data_start)

            for ind in range(datalen):
                update_progress(ind/datalen)
                row, col = [int(s) for s in lines[data_start[ind]+1].split()]
                tempdata = np.loadtxt(
                            lines[data_start[ind]+2:data_stop[ind]])

                which_hop = str(row) + ' ' + str(col)  

                total_hamdata[which_hop] = tempdata

        if return_data:
            return a, atom_pos, atom_names, atom_orbs, total_hamdata


















