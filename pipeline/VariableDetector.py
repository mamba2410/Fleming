from . import Utilities, Config, PeriodFinder

class VariableDetector:

    config = None
    n_sources = 0

    def __init__(self, config, n_sources, means, stds, medians, sorted_ids):
        self.config = config
        self.n_sources = n_sources
        self.sorted_ids = sorted_ids


    def get_variable_score(self, source_id, sorted_ids, sorted_stds):
        """
        Returns a score determining the variability of the star.

        Standard deviation search.
        
        Score is the standard deviation of the source compared
        to the median of the surrounding stars (in index-space)

        Parameters
        ----------

        index: int
            Index of variable star to check?

        adjusted: bool, optional
            Has the source light curve been adjusted?

        """
        
        source_index = np.where(self.sorted_ids == source_id)[0][0]
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim > self.n_sources:
            ulim = self.n_sources


        ## Median of standard deviations
        median_std = np.median(sorted_stds[llim:ulim])
        variable_score = sorted_stds[source_index]/median_std - 1

        return variable_score
        
        
    ## TODO: move signal to noise out of here
    def std_dev_search(self, adjusted=True):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable.

        Parameters
        ----------
        
        adjusted: bool, optional
            Has the source light curve been adjusted?


        Returns
        -------

        self.variable_ids: numpy array
            Numpy array of all ids which are deemed to be variable candidates.

        """

        if adjusted:
            sorted_stds = self.adjusted_source_stds[brightness_indices]
            sorted_medians = self.adjusted_source_medians[brightness_indices]
        else:
            sorted_stds = self.source_stds[brightness_indices]
            sorted_medians = self.source_medians[brightness_indices]


        self.variable_scores = np.zeros(self.n_sources)
        mins = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i, path, source_id in Utilities.loop_variables(self.config, self.source_ids, adjusted=adjusted):
            self.variable_scores[i] = self.get_variable_score(
                    source_id, sorted_ids, sorted_stds)
            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)
            mins[i] = np.min(lc['counts'])

        signal_to_noise = sorted_medians/sorted_stds
        self.variable_mask = np.where(
                (self.variable_scores > self.config.variability_threshold)
                & (self.variable_scores < self.config.variability_max)
                & (signal_to_noise > self.config.min_signal_to_noise)
                & (mins > 0)
                )[0]

        self.variable_ids_s = np.copy(self.source_ids[self.variable_mask])

            #else:
            #    print("[DEBUG] Rejecting, score of {}".format(self.variable_scores[i]))
                
                ## TODO: Check if catalogue is always ordered in ids
                #row_index = np.where(cat['id'].data==variable_id)

                ## TODO: Build table? Join with existing data
                #if 'RA' in cat.colnames:
                #    self.results_table.add_row([
                #        int(variable_id),
                #        cat['xcentroid'][row_index],
                #        cat['ycentroid'][row_index],
                #        variability,
                #        cat['RA'][row_index],
                #        cat['DEC'][row_index]])
                #else:
                #    self.results_table.add_row([
                #        int(variable_id),
                #        cat['xcentroid'][row_index],
                #        cat['ycentroid'][row_index], 
                #        variability,
                #        0,
                #        0])
        

        ## TODO: write results_table to file?

        print("[VariableDetector] std search: Found {} variables out of {} sources"
                .format(len(self.variable_ids), self.n_sources))

        return self.variable_ids_s


    def amplitude_search(self, source_ids=None, adjusted=True):
        """

        """

        if source_ids == None:
            source_ids = self.source_ids

        pf = PeriodFinder(self.config)

        amplitude_score = np.zeros(self.n_sources)

        for i, path, source_id in Utilities.loop_variables(self.config, source_ids, adjusted=adjusted):

            ## TODO: Make period range config variable
            _P, _P_err, A, A_err, _phi, _offset = pf.period_search(source_id, path)
            amplitude_score[i] = A/A_err

        variable_mask = np.where(amplitude_score > self.config.amplitude_score_threshold)[0]

        variable_ids = np.copy(source_ids[variable_mask])
        self.variable_ids_a = variable_ids

        print("[VariableDetector] amplitude search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def get_variables(self):
        """
        Function the user should call to get final IDs of sources deemed variable
        candiates.

        Returns
        -------

        variable_ids: numpy array
            IDs of each source considered variable

        """

        ## TODO: move this (to flux finder?)
        ## Remove any cosmic rays in the light curve
        for i, path, source_id in Utilities.loop_variables(self.config, self.source_ids, adjusted=True):
            curve = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()

            ## TODO: sigma clip
            if self.remove_cosmics(curve, stds[i]):
                #print("[DataAnalyser] Removed cosmics on id {}".format(source_id))
                np.savetxt(path, curve)

        ## Get variables judged by different methods
        variable_ids_s = self.std_dev_search(adjusted=True)
        variable_ids_a = self.amplitude_search(adjusted=True)

        ## Stars which apppear in both are returned
        variable_ids = np.intersect1d(variable_ids_s, variable_ids_a)


        return variable_ids


    def set_statistics(self, means, stds, medians, adjusted=False):
        if adjusted:
            self.source_means = means
            self.source_stds  = stds
            self.source_medians = medians
        else:
            self.adjusted_source_means = means
            self.adjusted_source_stds  = stds
            self.adjusted_source_medians = medians
