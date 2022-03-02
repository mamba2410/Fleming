import Utilities

import numpy as np
import os

import matplotlib.pyplot as plt

class PeriodFinder:

    def __init__(self, config):
        self.config = config


    ## Finds the global minimum
    def period_search(self, source_id, path,
            period_min = 1*3600, # 1 hour
            period_max = 5*3600, # 5 hours
            n_samples  = 2e3,    # 1000 samples
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

        omega_next = omegas[idx_min+1]
        chi2_next  = chi2[idx_min+1]

        chi2_orig = chi2
        omegas_orig = omegas



        ## TODO: Put this in a loop until some tolerance?
        iteration = 1
        while (np.max(chi2) - chi2_min) > self.config.period_chi2_range \
              and iteration < self.config.period_max_iterations:

            print("[PeriodFinder] Iteration: {:02}".format(iteration))

            ## Approximate chi2 landscape as linear, for rough domain searches
            approx_gradient = np.sqrt(chi2_next - chi2_min)/(omega_next - omega_min)

            ## Estimation of step we should take to get \deltachi2 ~1
            approx_halfwidth = 1/approx_gradient

            omegas = np.linspace(
                    omega_min - self.config.period_width_adjustment*approx_halfwidth,
                    omega_min + self.config.period_width_adjustment*approx_halfwidth,
                    int(n_samples)
                )

            ## Second approximation, much closer to the minimum
            ## Need to get close enough to be able to approximate to a parabola
            _params, chi2 = self.search_omegas(counts, time, errors, omegas)
            idx_min = np.argmin(chi2)
            omega_min = omegas[idx_min]
            chi2_min  = chi2[idx_min]

            plt.plot(2*np.pi/omegas, chi2)
            fname = "chi2_i{:02}_{}_{}{:04}{}".format(
                iteration,
                self.config.image_prefix,
                self.config.identifier,
                source_id,
                self.config.plot_file_extension
            )
            plt.title("Iteration {:02}".format(iteration))
            plt.savefig(os.path.join(self.config.periods_dir, fname))
            plt.close()

            iteration += 1

        ## Assume we are good enough
        ## TODO: End of loop

        idx_dchi2 = np.where(chi2[idx_min:] > chi2_min + 1)[0][0]
        idx_dchi2 += idx_min
        
        c = (chi2[idx_dchi2] - chi2_min)/(omegas[idx_dchi2] - omega_min)**2
        omega_min_err = np.sqrt(1/c)


        period_min = 2*np.pi/omega_min
        period_min_err = (omega_min_err/omega_min) * period_min

        ## Sine curve of found period
        B, C, S = _params[idx_min]
        attempted_fit = B + C*np.cos(omega_min*time) + S*np.sin(omega_min*time)

        plt.scatter(time, counts)
        plt.plot(time, attempted_fit, color="red")
        fname = "compare_{}_{}{:04}_P{:6g}{}".format(
                self.config.image_prefix,
                self.config.identifier,
                source_id,
                period_min,
                self.config.plot_file_extension
            )
        plt.title("Comparison of found period and light curve for id{:04}".format(source_id))
        plt.xlabel("Time [s]")
        plt.ylabel("Normalised brightness [arb. u.]")
        plt.savefig(os.path.join(self.config.periods_dir, fname))
        plt.close()


        params, _H, H_inv = self.periodogram_fit_func(counts, time, errors, omega_min)
        B, C, S = params
        C_err = H_inv[1,1]
        S_err = H_inv[2,2]
        A = np.sqrt(C**2 + S**2)

        amplitude_err = np.sqrt( ((C/A)*C_err)**2 + ((S/A)*S_err)**2 )
        amplitude = A

        ## Phase shift for a sin wave
        phi = np.pi/2-np.arctan(S/C)

        ## Plot initial chi2 distribution
        plt.plot(2*np.pi/omegas_orig, chi2_orig)
        fname = "chi2_{}_{}{:04}_P{:6g}{}".format(
                self.config.image_prefix,
                self.config.identifier,
                source_id,
                period_min,
                self.config.plot_file_extension
            )
        plt.title("$\chi^2$ plot for source {:04}".format(source_id))
        plt.xlabel("Period [s]")
        plt.ylabel("$\chi^2$ [arb. u.]")
        plt.savefig(os.path.join(self.config.periods_dir, fname))
        plt.close()


        return period_min, period_min_err, amplitude, amplitude_err, phi, B



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


