import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

h2eV = 27.2114
efermi = 2.59236423E-01*h2eV

nband = 20

eigenvalues = nc.Dataset('runo_DS1_EIG.nc')['Eigenvalues'][0,:,:]*h2eV-efermi
nkpt = eigenvalues.shape[0]
nband = eigenvalues.shape[1]

fig,ax=plt.subplots(1,2,subplot_kw={'sharex':True,'sharey':True},
        gridspec_kw={'wspace':0.1,'hspace':0.3})
fig.set_size_inches(10,5,forward=True)

for i in range(nband):
    ax[0].plot(eigenvalues[:,i],ls='',marker='o',ms=3,mew=1,mfc='r',mec='k')
ax[0].plot([0,nkpt],[0,0],'k',ls=':')

img = mpimg.imread('ref_ebands.png')
ax[1].imshow(img)

ax[0].axis([-1,nkpt,-15,10])
ax[1].set_yticks([])

plt.savefig('comparison_0.05.pdf',format='pdf',dpi=150)
plt.show()
