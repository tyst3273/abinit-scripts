import numpy as np

def write_strufac(params,data):
    
    sfac = data.structure_factors
    freq = np.zeros((sfac.shape[0],sfac.shape[1]))

    for q in range(data.num_Qpoints):
        freq[:,q] = data.frequencies[f'{q}']
    
    freq = np.append(freq,sfac,axis=0)
    np.savetxt('structure_factors',freq,header='freq are in meV\nfreq(b=1,q=1) freq(b=1,q=2) freq(b=1,q=3) ...\nfreq(b=2,q=1) freq(b=2,q=2) freq(b=2,q=3) ...\n...\nsfac(b=1,q=1) sfac(b=1,q=2) sfac(b=1,q=3) ...\nsfac(b=2,q=1) sfac(b=2,q=2) sfac(b=2,q=3) ...\n...')
