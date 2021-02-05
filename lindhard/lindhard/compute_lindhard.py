import numpy as np
import matplotlib.pyplot as plt
import h5py 


### some parameters
t_smear = 300 # kelvin
k_prec = 6

### read in data
in_f_name = 'eigs.hdf5'
in_file = h5py.File(in_f_name,'r')
eigenvals = in_file['eigenvals'][:,:]
kpts = in_file['kpts'][:,:]
in_file.close()
n_kpts = kpts.shape[0]
n_bands = eigenvals.shape[1]

### k-vector step sizes
dh = np.unique(kpts[:,0])
dh = dh[np.argsort(dh)]
dh = np.round(dh[1]-dh[0],k_prec)
dk = np.unique(kpts[:,0])
dk = dk[np.argsort(dk)]
dk = np.round(dk[1]-dk[0],k_prec)
dl = np.unique(kpts[:,0])
dl = dl[np.argsort(dl)]
dl = np.round(dl[1]-dl[0],k_prec)

k_str = []
for k in range(n_kpts):
    print(k)
