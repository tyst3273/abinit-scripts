import netCDF4 as nc
import h5py


n_dirs = 100


h2eV = 27.2114

n_buff = 4
e_fermi = 2.59236423e-01 # hartree, should be the same for all files

out_f_name = 'eigs.hdf5'
out_file = h5py.File(out_f_name,'w')

for i in range(1,n_dirs+1):

    print(i)

    eig_dir = f'eigs_{i}'
    in_f_name = eig_dir+'/'+'eigso_DS1_EIG.nc'
    dat = nc.Dataset(in_f_name)

    if i == 1:
        n_bands = dat['Eigenvalues'][0,0,:].shape[0]-n_buff
        n_kpts = dat['Eigenvalues'][0,:,0].shape[0]
        eigenvals = out_file.create_dataset('eigenvals',
                [n_kpts*n_dirs,n_bands])
        kpts = out_file.create_dataset('kpts',[n_kpts*n_dirs,3])

    eigenvals[(i-1)*n_kpts:i*n_kpts,:] = dat['Eigenvalues'][
            0,:,:n_bands]-e_fermi 
    kpts[(i-1)*n_kpts:i*n_kpts,:] = dat['Kptns'][:,:]

eigenvals[:,:] = eigenvals[:,:]*h2eV

    

out_file.close()

