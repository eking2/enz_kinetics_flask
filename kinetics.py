import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from io import StringIO

plt.rcParams['font.sans-serif'] = 'Arial'

def make_assay_df(string):

    '''convert space delimited input to df'''

    # expect data do have concs in column 1
    # trials in next set of columns

    # read in string to df, no headers
    inp = StringIO(string)
    df = pd.read_table(inp, delim_whitespace=True, header=None)

    # set columns names
    col_names = [f't{i+1}' for i in range(len(df.columns) - 1)]
    col_names = ['cofa_conc_mM'] + col_names
    df.columns = col_names

    # convert units
    # from mA/min to A/s
    df.iloc[:, 1:] = df.iloc[:, 1:] * 1/60 * 1/1000

    # convert to tidy
    # columns currently [substrate concs, t1, t2...]
    # melt to [trial, substrate concs, slope]
    df = pd.melt(frame=df, id_vars = 'cofa_conc_mM', var_name='trial', value_name='slope_a_s')

    return df


def make_enz_dict(rxn_vol, enz_vol, enz_conc, mol_wt, dil, ext):
    
    '''convert units for michaelis menten calc'''
    
    # mass of enzyme in reaction
    # mg/ml (same as g/L) * 1/df * ul/ul = g/L
    enz_rxn_conc_g_l = enz_conc * 1/dil * enz_vol/rxn_vol

    # g/L * mol/g * 10^6umol/mol
    enz_rxn_uM = enz_rxn_conc_g_l * 1/mol_wt * 10**6
    enz_rxn_mM = float(enz_rxn_uM / 1000)

    # pathlength from rxn vol
    # 200uL -> 0.5cm, 100uL -> 0.25cm
    pathlen = int(rxn_vol) / 400

    enz_dict = {'ext' : float(ext),         # mM^-1 cm^-1
                'pathlen' : pathlen,        # cm
                'enz_rxn_mM' : enz_rxn_mM}  # mM

    return enz_dict


def velocity(sub, enz, kcat, km):

    '''michaelis menten equation

    sub : float
        concentration of substrate in mM
    enz : float
        concentration of enzyme in mM
    '''

    return (kcat * enz * sub) / (km + sub)


class kinetics_calc:

    def __init__(self, assay_df, ext, pathlen, enz_rxn_mM):

        self.assay_df = assay_df
        self.ext = ext
        self.pathlen = pathlen
        self.enz_rxn_mM = enz_rxn_mM

        self.abs_to_vel()


    def abs_to_vel(self):

        '''abs/sec to v0/sec (mM/sec)'''

        # A/sec * mM
        # cm cancels
        self.assay_df['v0_mM_s'] = self.assay_df['slope_a_s'] / (self.ext * self.pathlen)


    def fit_mm(self):

        '''fit data to michaelis menten'''

        # fix the enzyme conc with lambda
        popt, pcov = curve_fit(lambda sub, kcat, km: velocity(sub, self.enz_rxn_mM, kcat, km), 
                self.assay_df['cofa_conc_mM'], self.assay_df['v0_mM_s'])

        perr = np.sqrt(np.diag(pcov))

        # coefficients and errors (kcat, km)
        return popt, perr


    def get_rsq(self, popt):

        '''get rsq for data to fit curve'''
        
        residuals = self.assay_df['v0_mM_s'] - velocity(self.assay_df['cofa_conc_mM'], self.enz_rxn_mM, popt[0], popt[1])
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((self.assay_df['v0_mM_s'] - np.mean(self.assay_df['v0_mM_s']))**2)
        r_sq = 1 - (ss_res / ss_tot)

        return r_sq


    def plot_mm(self, title=None):

        '''plot michaelis menten'''

        # get kinetics
        popt, perr = self.fit_mm()
        r_sq = self.get_rsq(popt)

        kcat, km = popt[0], popt[1]
        kcat_err, km_err = perr[0], perr[1]

        # range from lowest cofa conc to highest
        min_cofa = self.assay_df['cofa_conc_mM'].min()
        max_cofa = self.assay_df['cofa_conc_mM'].max()
        x = np.linspace(min_cofa, max_cofa, 1000)

        # best fit curve
        plt.plot(x, velocity(x, self.enz_rxn_mM, kcat, km) / self.enz_rxn_mM, color='C0', lw=2, label='Michaelis Menten')

        # km line
        plt.axvline(km, color='orange', ls='--', lw=1.5, alpha=0.6, label='$K_m$')

        # scatter experimental data 
        plt.scatter(self.assay_df['cofa_conc_mM'], self.assay_df['v0_mM_s'] / self.enz_rxn_mM, color='limegreen', edgecolor='k',
                    zorder=10, s=30, alpha=0.8)

        # annotate
        annotation = r'$k_{{cat}}$ = {:.3f} $\pm$ {:.2f} s$^{{-1}}$'.format(kcat, kcat_err)
        annotation += '\n'
        annotation += r'$K_m$ = {:.3f} $\pm$ {:.2f} mM'.format(km, km_err)
        annotation += '\n'
        annotation += r'$R^2$ = {:.3f}'.format(r_sq)

        plt.text(0.47, 0.13, annotation, transform=plt.gca().transAxes, size=14, linespacing=1.6)

        # clean up
        if title:
            plt.title(title, size=18)

        plt.ylabel('Rate (s$^{-1}$)', size=15)
        plt.xlabel('Substrate (mM)', size=15)
        plt.grid(alpha=0.2)
        plt.gca().set_axisbelow(True)
        plt.savefig('./static/plot.png', bbox_inches='tight', dpi=600)






