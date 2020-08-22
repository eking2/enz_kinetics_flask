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


def lin_velocity(sub, enz, cat_eff):

    '''linear michaelis menten when K_m >> S, cat_eff = kcat/km'''

    return enz * cat_eff * sub


class kinetics_calc:

    def __init__(self, assay_df, fit, ext, pathlen, enz_rxn_mM):

        self.assay_df = assay_df
        self.ext = ext
        self.pathlen = pathlen
        self.enz_rxn_mM = enz_rxn_mM
        self.fit = fit

        self.r_sq = None
        self.kcat = None
        self.kcat_err = None
        self.km = None
        self.km_err = None
        self.cat_eff = None
        self.cat_eff_err = None

        self.abs_to_vel()
        self.fit_mm()
        self.get_rsq()


    def abs_to_vel(self):

        '''abs/sec to v0/sec (mM/sec)'''

        # A/sec * mM
        # cm cancels
        self.assay_df['v0_mM_s'] = self.assay_df['slope_a_s'] / (self.ext * self.pathlen)


    def fit_mm(self):

        '''fit data to michaelis menten'''

        # coefficients and errors
        # (kcat, km) for hyperbolic
        # kcat/km for linear

        if self.fit == 'hyperbolic':
            # fix the enzyme conc with lambda
            popt, pcov = curve_fit(lambda sub, kcat, km: velocity(sub, self.enz_rxn_mM, kcat, km),
                    self.assay_df['cofa_conc_mM'], self.assay_df['v0_mM_s'])

            perr = np.sqrt(np.diag(pcov))

            self.kcat, self.km = popt[0], popt[1]
            self.kcat_err, self.km_err = perr[0], perr[1]

            # propagate error
            covar = pcov[0][-1]
            self.cat_eff = self.kcat / self.km
            self.cat_eff_err = np.absolute(self.cat_eff) * ( (self.kcat_err/self.kcat)**2 + (self.km_err/self.km)**2 + 2*covar/(self.kcat * self.km))

        else:
            popt, pcov = curve_fit(lambda sub, cat_eff: lin_velocity(sub, self.enz_rxn_mM, cat_eff),
                    self.assay_df['cofa_conc_mM'], self.assay_df['v0_mM_s'])

            perr = np.sqrt(np.diag(pcov))

            self.cat_eff = popt[0]
            self.cat_eff_err = perr[0]


    def get_rsq(self):

        '''get rsq for data to fit curve'''

        if self.fit == 'hyperbolic':
            residuals = self.assay_df['v0_mM_s'] - velocity(self.assay_df['cofa_conc_mM'], self.enz_rxn_mM, self.kcat, self.km)

        else:
            residuals = self.assay_df['v0_mM_s'] - lin_velocity(self.assay_df['cofa_conc_mM'], self.enz_rxn_mM, self.cat_eff)

        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((self.assay_df['v0_mM_s'] - np.mean(self.assay_df['v0_mM_s']))**2)
        self.r_sq = 1 - (ss_res / ss_tot)


    def plot_mm(self, title=None):

        '''plot michaelis menten'''

        # range from lowest cofa conc to highest
        min_cofa = self.assay_df['cofa_conc_mM'].min()
        max_cofa = self.assay_df['cofa_conc_mM'].max()
        x = np.linspace(min_cofa, max_cofa, 1000)

        if self.fit == 'hyperbolic':

            # best fit curve
            # plot rate on y instead of velocity, divide by enz conc
            plt.plot(x, velocity(x, self.enz_rxn_mM, self.kcat, self.km) / self.enz_rxn_mM, color='C0', lw=2, label='Michaelis Menten')

            # km line
            plt.axvline(self.km, color='orange', ls='--', lw=1.5, alpha=0.6, label='$K_m$')

            # annotate
            annotation = r'$k_{{cat}}$ = {:.3f} $\pm$ {:.2f} s$^{{-1}}$'.format(self.kcat, self.kcat_err)
            annotation += '\n'
            annotation += r'$K_M$ = {:.3f} $\pm$ {:.2f} mM'.format(self.km, self.km_err)
            annotation += '\n'

            # text too long to fit on plot
            #annotation += r'$k_{{cat}}$ / $K_m$ = {:.3f} $\pm$ {:.2f} s$^{{-1}}$ mM$^{{-1}}$'.format(self.cat_eff, self.cat_eff_err)
            #annotation += '\n'
            annotation += r'$R^2$ = {:.3f}'.format(self.r_sq)

            plt.text(0.47, 0.13, annotation, transform=plt.gca().transAxes, size=14, linespacing=1.6)

        else:

            # plot linear fit
            plt.plot(x, lin_velocity(x, self.enz_rxn_mM, self.cat_eff) / self.enz_rxn_mM, color='C0', lw=2, label='Linear Michaelis Menten')

            annotation = r'$k_{{cat}} / K_M$ = {:.2f} $\pm$ {:.1f} s$^{{-1}}$ mM$^{{-1}}$'.format(self.cat_eff, self.cat_eff_err)
            annotation += '\n'
            annotation += r'$R^2$ = {:.3f}'.format(self.r_sq)

            plt.text(0.37, 0.13, annotation, transform=plt.gca().transAxes, size=12, linespacing=1.6)


        # scatter experimental data
        plt.scatter(self.assay_df['cofa_conc_mM'], self.assay_df['v0_mM_s'] / self.enz_rxn_mM, color='limegreen', edgecolor='k',
                    zorder=10, s=30, alpha=0.8)

        # clean up
        if title:
            plt.title(title, size=18)

        plt.ylabel('Rate (s$^{-1}$)', size=15)
        plt.xlabel('Substrate (mM)', size=15)
        plt.grid(alpha=0.2)
        plt.gca().set_axisbelow(True)
        plt.savefig('./static/plot.png', bbox_inches='tight', dpi=600)

        plt.close()


    def save_output(self):

        '''output dict with solved kinetic parameters, save dataframe'''

        output = {}

        # convert to regular float so yaml is human readable
        # otherwise prints extra
        output['kcat'] = float(self.kcat) if self.kcat else None
        output['kcat_err'] = float(self.kcat_err) if self.kcat_err else None
        output['km'] = float(self.km) if self.km  else None
        output['km_err'] = float(self.km_err) if self.km_err else None
        output['cat_eff'] = float(self.cat_eff)
        output['cat_eff_err'] = float(self.cat_eff_err)
        output['r_sq'] = float(self.r_sq)
        output['fit'] = self.fit

        self.assay_df[['trial', 'cofa_conc_mM', 'slope_a_s', 'v0_mM_s']].to_csv('static/assay_df.csv', index=False)

        return output




