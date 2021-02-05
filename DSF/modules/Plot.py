import matplotlib.pyplot as plt
import numpy as np

class plot_dispersions:
    
    def __init__(self,params,data):

        fig, ax = plt.subplots()
        fig.set_size_inches(6,4,forward=True)
        fig.tight_layout(pad=4)
        q_array = np.arange(params.num_Qpoints)
        for q in range(params.num_Qpoints):
            for b in range(data.num_bands):
                ax.plot(q_array[q],data.frequencies[f'{q}'][b],
                        marker='o',ms=4,mfc='none',mec='k',mew=2)

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='medium')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)
        ax.set_ylabel('E (meV)',
                        labelpad=5.0,
                        fontweight='normal',
                        fontsize='large')
        ax.set_xlabel('Q index',
                        labelpad=8.0,
                        fontweight='normal',
                        fontsize='large')

        fig.suptitle('Phonon Dispersion',y=0.98,fontsize='large')
        if params.save_figs == True:
            plt.savefig('phonon-dispersion.png',fmt='png',dpi=300,transparent=False)
        plt.show()


class plot_structure_factors:

    def __init__(self,params,data):

        fig, ax = plt.subplots()
        fig.set_size_inches(6,4,forward=True)
        fig.tight_layout(pad=4)

        sfac = data.structure_factors
        q_array = np.arange(params.num_Qpoints)

        if params.sum_degenerate_bands == True:

            print('\n\tSumming degenerate bands before plotting\n')
            for q in range(params.num_Qpoints):
                energy = np.unique(np.round(data.frequencies[f'{q}'],params.degeneracy_tolerance_decimals))
                summed_sfac = np.zeros(len(energy))

                for i in range(len(energy)):
                    summed_sfac[i] = sfac[np.argwhere(
                        np.round(data.frequencies[f'{q}'],params.degeneracy_tolerance_decimals) == energy[i]),q].sum()

                for b in range(len(energy)):
                    ax.plot(q_array[q],energy[b],
                            marker='o',ms=summed_sfac[b]*5,mfc='none',mec='b',mew=1)

        else:

            for q in range(params.num_Qpoints):
                for b in range(data.num_bands):
                    ax.plot(q_array[q],data.frequencies[f'{q}'][b],
                            marker='o',ms=sfac[b,q]*5,mfc='none',mec='b',mew=1)

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='medium')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)
        ax.set_ylabel('E (meV)',
                        labelpad=5.0,
                        fontweight='normal',
                        fontsize='large')
        ax.set_xlabel('Q index',
                        labelpad=8.0,
                        fontweight='normal',
                        fontsize='large')
        fig.suptitle('Structure Factors',y=0.98,fontsize='large')
        if params.save_figs == True:
            plt.savefig('phonon-structure-factors.png',fmt='png',dpi=300,transparent=False)
        plt.show()




