import Utilities
import numpy as np

import matplotlib.pyplot as plt

class PeriodFinder:

    def __init__(self, config):
        self.config = config


    def period_search(self, source_id, path,
            period_min = 1*3600, # 1 hour
            period_max = 5*3600, # 5 hours
            n_samples  = 1e4,    # 10_000 samples
            adjusted=True):
        
        curve = np.genfromtxt(path, dtype = self.config.light_curve_dtype).transpose()
        time = curve['time']
        counts = curve['counts']
        errors = curve['counts_err']
        
        periods = np.linspace(period_min, period_max, int(n_samples))
        omegas = 2*np.pi/periods


        ## Initial search

        _params, chi2 = self.search_omegas(counts, time, errors, omegas)
        idx_min = np.argmin(chi2)
        omega_min = omegas[idx_min]
        chi2_min  = chi2[idx_min]

        plt.plot(2*np.pi/omegas, chi2)
        plt.savefig(self.config.output_dir + "/TEST1.jpg")
        plt.close()

        omega_next = omegas[idx_min+1]
        chi2_next  = chi2[idx_min+1]

        ## TODO: Put this in a loop until some tolerance?

        ## Approximate chi2 landscape as linear, for rough domain searches
        approx_gradient = (chi2_next - chi2_min)/(omega_next - omega_min)

        ## Estimation of step we should take to get \deltachi2 ~1
        approx_halfwidth = 1/approx_gradient

        ## How much should we scale this width to capture the whole minimum?
        width_adjustment = 2

        omegas = np.linspace(
                omega_min - width_adjustment*approx_halfwidth,
                omega_min + width_adjustment*approx_halfwidth,
                int(n_samples)
            )

        ## Second approximation, much closer to the minimum
        _params, chi2 = self.search_omegas(counts, time, errors, omegas)
        idx_min = np.argmin(chi2)
        omega_min = omegas[idx_min]
        chi2_min  = chi2[idx_min]

        ## Assume we are good enough

        idx_dchi2 = np.where(chi2[idx_min:] > chi2_min + 1)[0][0]
        idx_dchi2 += idx_min
        
        c = (chi2[idx_dchi2] - chi2_min)/(omegas[idx_dchi2] - omega_min)**2
        omega_min_err = np.sqrt(1/c)

        plt.plot(2*np.pi/omegas, chi2)
        plt.savefig(self.config.output_dir + "/TEST2.jpg")
        plt.close()


        period_min = 2*np.pi/omega_min
        period_min_err = (omega_min_err/omega_min) * period_min

        return period_min, period_min_err



    def search_omegas(self, counts, time, errors, omegas):
        n_samples = len(omegas)
        chi2 = np.zeros(n_samples)
        params = np.zeros((n_samples, 3))

        for i,omega in enumerate(omegas):
            params[i], _H, _H_inv = self.periodogram_fit_func(counts, time, errors, omega)
            chi2[i] = self.periodogram_chi2_(counts, time, errors, omega, params[i])

        return params, chi2
            
        


    def periodogram_fit_func(self, data, time, error, omega):
        s2r = 1/error**2
        s = np.sin(omega*time)
        c = np.cos(omega*time)
        
        H = np.array([
            [np.sum(s2r), np.sum(c*s2r), np.sum(s*s2r)],
            [np.sum(c*s2r), np.sum(c*c*s2r), np.sum(c*s*s2r)],
            [np.sum(s*s2r), np.sum(c*s*s2r), np.sum(s*s*s2r)]
        ])
        
        corr = np.array([
            [np.sum(data*s2r)],
            [np.sum(data*c*s2r)],
            [np.sum(data*s*s2r)]
        ])
        
        H_inv = np.linalg.inv(H)
        fit_params = np.reshape(np.matmul(H_inv, corr), 3)
        return fit_params, H, H_inv


    def periodogram_chi2(self, data, time, error, omega, B, C, S):
        F = lambda t: B + C*np.cos(omega*t) + S*np.sin(omega*t)
        chi2 = ((data-F(time))/error)**2
        chi2 = np.sum(chi2)
        return chi2
    
    
    def periodogram_chi2_(self, data, time, error, omega, params):
        return self.periodogram_chi2(data, time, error, omega, params[0], params[1], params[2])


