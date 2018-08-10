#! /usr/bin/env python
'''
generate initial guess of u matrix from small cluster.

'''

import numpy as np
import h5py
import sys
from hfb_hubbard_2 import utils

DEBUG = True


def replicate_u_matrix( u_imp, Nbasis, Nimp, mtype = np.float64):

    #Number of copies of impurity cluster in total lattice
    Ncopies = Nbasis / Nimp

    #Copy impurity u-matrix across entire lattice
    u_mat_replicate = np.zeros( [2*Nbasis,2*Nbasis],dtype=mtype)
    for cpy in xrange(Ncopies):
        start1 = cpy * Nimp
        end1 = start1 + Nimp
        start2 = start1 + Nbasis
        end2 = end1 + Nbasis

        u_mat_replicate[start1:end1, start1:end1] = u_imp[:Nimp, :Nimp]
        u_mat_replicate[start2:end2, start1:end1] = u_imp[Nimp:Nimp+Nimp, :Nimp]
        u_mat_replicate[start1:end1, start2:end2] = u_imp[:Nimp, Nimp:Nimp+Nimp]
        u_mat_replicate[start2:end2, start2:end2] = u_imp[Nimp:Nimp+Nimp, Nimp:Nimp+Nimp]
    return u_mat_replicate

def array2matrix( array, Nimp, mtype  = np.float64):

    mat = np.zeros([Nimp*2, Nimp*2], mtype)
    
    mat[np.triu_indices(Nimp * 2)] = array
    mat = mat + mat.conj().T
    mat[np.diag_indices(Nimp * 2)] *= 0.5
    
    return mat

def matrix2array( mat, Nimp):
        
    array = mat[np.triu_indices(Nimp * 2)]
    return array

np.set_printoptions(8, linewidth =1000)

f1_name = sys.argv[1]
f1 = h5py.File(f1_name)

# global mu
mu_1 = f1['globalMu.r'].value.reshape(f1['globalMu.s'].value)
# umatrix
umat_1 = f1['umatrix.r'].value.reshape(f1['umatrix.s'].value)
# rdm
rdm_1 = f1['rdm1.r'].value.reshape(f1['rdm1.s'].value)
# mu for correlation calc
muc_1 = f1['muc.r'].value.reshape(f1['muc.s'].value)

if DEBUG:
    print f1.keys()
    print 
    print "mu_1"
    print mu_1
    print
    print "umat_1"
    print umat_1
    print
    print "muc_1"
    print muc_1

np.set_printoptions(3,linewidth= 2000)
R = rotmat_1 = f1['rotationmatrix.r'].value.reshape(f1['rotationmatrix.s'].value)
umat = umat_1

lbasis, lactorb = R.shape
Nbasis = lbasis / 2
Nimp = lactorb / 4
limp = 2 * Nimp
print "umat"
print umat


Ix = 2
Iy = 2
I = [Ix, Iy]
U = 4.0
Filling = 0.4

'''
um = np.zeros([2*Ix*Iy,2*Ix*Iy])
if (True):
    um_u, um_d, um_o = utils.AFInitGuess(I, U, Filling, polar = None, bogoliubov = True, rand = 0.001)
    #um_u, um_d, um_o = PMInitGuess(I, U, Filling, bogoliubov = True, rand = 0.001)
    um[:Ix*Iy, :Ix*Iy] = um_u
    um[Ix*Iy:, Ix*Iy:] = -um_d
    um[:Ix*Iy, Ix*Iy:] = um_o
    um[Ix*Iy:, :Ix*Iy] = um_o.conjugate().T
print "um"
print um
'''

u22 = umat
u42_uu = np.zeros((8,8))
u42_uu[np.ix_([0,1,4,5], [0,1,4,5])] = u22[:4, :4]
u42_uu[np.ix_([2,3,6,7], [2,3,6,7])] = u22[:4, :4]
#u42_uu[np.ix_([1,2,5,6], [1,2,5,6])] = u22[:4, :4]
#u42_uu[np.ix_([3,0,7,4], [3,0,7,4])] = u22[:4, :4]

u42_ud = np.zeros((8,8))
u42_ud[np.ix_([0,1,4,5], [0,1,4,5])] = u22[:4, 4:]
u42_ud[np.ix_([2,3,6,7], [2,3,6,7])] = u22[:4, 4:]
#u42_ud[np.ix_([1,2,5,6], [1,2,5,6])] = u22[:4, 4:]
#u42_ud[np.ix_([3,0,7,4], [3,0,7,4])] = u22[:4, 4:]

u42_dd = np.zeros((8,8))
u42_dd[np.ix_([0,1,4,5], [0,1,4,5])] = u22[4:, 4:]
u42_dd[np.ix_([2,3,6,7], [2,3,6,7])] = u22[4:, 4:]
#u42_dd[np.ix_([1,2,5,6], [1,2,5,6])] = u22[4:, 4:]
#u42_dd[np.ix_([3,0,7,4], [3,0,7,4])] = u22[4:, 4:]

u42_imp = np.zeros((16,16))
u42_imp[:8, :8] = u42_uu
u42_imp[:8, 8:] = u42_ud
u42_imp[8:, 8:] = u42_dd
u42_imp[8:, :8] = u42_ud.conj().T

print "u42_imp uu"
print u42_imp[:8, :8]
print "u42_imp dd"
print u42_imp[8:, 8:]

u42_imp_new = (u42_imp + u42_imp.conj().T) * 0.5
np.save("um0_4x2.npy", u42_imp_new)
#print np.linalg.norm(u42_imp_new - u42_imp)

#u22 = np.arange(1,17).reshape(4,4)
#u42 = np.zeros((8,8))
#u42[np.ix_([0,1,4,5], [0,1,4,5])] = u22

#print u22
#print u42

exit()


'''
print "umat"
print umat
umat_full = replicate_u_matrix( umat, Nbasis, Nimp)
umat_add = np.dot(R.conj().T, np.dot(umat_full, R))
print "umat act"
print umat_add
'''
