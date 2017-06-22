import warnings
import pickle
import math

import numpy as np
import scipy.linalg as linalg
import scipy.ndimage as ndi
import scipy.ndimage.interpolation as sinterp

import pyklip.covars as covars

# emcee more MCMC sampling
import emcee

# plotting tools
import matplotlib
import matplotlib.pylab as plt
import corner



class FMAstrometry(object):
    """
    Base class to perform astrometry on direct imaging data_stamp using MCMC and GP regression
    """
    def __init__(self, guess_sep, guess_pa, fitboxsize):
        """
        Initilaizes the FMAstrometry class

        Args:
            guess_sep: the guessed separation (pixels)
            guess_pa: the guessed position angle (degrees)
            fitboxsize: fitting box side length (pixels)

        Returns:

        """
        # store initailization
        self.guess_sep = guess_sep
        self.guess_pa = guess_pa
        self.fitboxsize = fitboxsize

        # derive delta RA and delta Dec
        # in pixels
        self.guess_RA_offset = self.guess_sep * np.sin(np.radians(self.guess_pa))
        self.guess_Dec_offset = self.guess_sep * np.cos(np.radians(self.guess_pa))

        # stuff that isn't generated yet
        # stamps of the data_stamp and the forward model
        self.fm_stamp = None # Forward Model
        self.padding = 0 # padding for FM. You kinda need this to shift the FM around
        self.data_stamp = None # Data
        self.noise_map = None # same shape as self.data_stamp
        self.data_stamp_RA_offset = None # RA offset of data_stamp (in pixels)
        self.data_stamp_Dec_offset = None # Dec offset (in pixels)
        self.data_stamp_RA_offset_center = None # RA offset of center pixel (stampsize // 2)
        self.data_stamp_Dec_offset_center = None # Dec offset of center pixel (stampsize // 2)

        # guess flux (a hyperparameter)
        self.guess_flux = None

        # covariance paramters. Use the covariance initilizer function to initilize them
        self.covar = None
        self.covar_param_guesses = None
        self.covar_param_labels = None
        self.include_readnoise = False

        # MCMC fit params
        self.bounds = None
        self.sampler = None

        # best fit
        self.RA_offset = None
        self.RA_offset_1sigma = None
        self.Dec_offset = None
        self.Dec_offset_1sigma = None
        self.flux = None
        self.flux_1sigma = None
        self.covar_param_bestfits = None
        self.covar_param_1sigma = None

    def generate_fm_stamp(self, fm_image, fm_center=None, fm_wcs=None, extract=True, padding=5):
        """
        Generates a stamp of the forward model and stores it in self.fm_stamp
        Args:
            fm_image: full imgae containing the fm_stamp
            fm_center: [x,y] center of image (assuing fm_stamp is located at sep/PA) corresponding to guess_sep and guess_pa
            fm_wcs: if not None, specifies the sky angles in the image. If None, assume image is North up East left
            extract: if True, need to extract the forward model from the image. Otherwise, assume the fm_stamp is already
                    centered in the frame (fm_image.shape // 2)
            padding: number of pixels on each side in addition to the fitboxsize to extract to pad the fm_stamp
                        (should be >= 1)

        Returns:

        """
        # cheeck the padding to make sure it's valid
        if not isinstance(padding, int):
            raise TypeError("padding must be an integer")
        if padding < 1:
            warnings.warn("Padding really should be >= 1 pixel so we can shift the FM around", RuntimeWarning)
        self.padding = padding


        if extract:
            if fm_wcs is not None:
                raise NotImplementedError("Have not implemented rotation using WCS")

            # image is now rotated North up east left
            # find the location of the FM
            thistheta = np.radians(self.guess_pa + 90)
            psf_xpos = self.guess_sep * np.cos(thistheta) + fm_center[0]
            psf_ypos = self.guess_sep * np.sin(thistheta) + fm_center[1]

        else:
            # PSf is already cenetered
            psf_xpos = fm_image.shape[1]//2
            psf_ypos = fm_image.shape[0]//2

        # now we found the FM in the image, extract out a centered stamp of it
        # grab the coordinates of the image
        stampsize = 2 * self.padding + self.fitboxsize # full stamp needs padding around all sides
        x_stamp, y_stamp = np.meshgrid(np.arange(stampsize * 1.) - stampsize //2,
                                       np.arange(stampsize * 1.) - stampsize// 2)

        x_stamp += psf_xpos
        y_stamp += psf_ypos

        # zero nans because it messes with interpolation
        fm_image[np.where(np.isnan(fm_image))] = 0

        fm_stamp = ndi.map_coordinates(fm_image, [y_stamp, x_stamp])
        self.fm_stamp = fm_stamp



    def generate_data_stamp(self, data, data_center, data_wcs=None, noise_map=None, dr=4, exclusion_radius=10):
        """
        Generate a stamp of the data_stamp ~centered on planet and also corresponding noise map
        Args:
            data: the final collapsed data_stamp (2-D)
            data_center: location of star in the data_stamp
            data_wcs: sky angles WCS object. To rotate the image properly [NOT YET IMPLMETNED]
                      if None, data_stamp is already rotated North up East left
            noise_map: if not None, noise map for each pixel in the data_stamp (2-D).
                        if None, one will be generated assuming azimuthal noise using an annulus widthh of dr
            dr: width of annulus in pixels from which the noise map will be generated
            exclusion_radius: radius around the guess planet location which doens't get factored into noise estimate

        Returns:

        """
        # rotate image North up east left if necessary
        if data_wcs is not None:
            # rotate
            raise NotImplementedError("Rotating based on WCS is not currently implemented yet")

        xguess = -self.guess_RA_offset + data_center[0]
        yguess = self.guess_Dec_offset + data_center[1]

        # round to nearest pixel
        xguess_round = int(np.round(xguess))
        yguess_round = int(np.round(yguess))

        # get index bounds for grabbing pixels from data_stamp
        ymin = yguess_round - self.fitboxsize//2
        xmin = xguess_round - self.fitboxsize//2
        ymax = yguess_round + self.fitboxsize//2 + 1
        xmax = xguess_round + self.fitboxsize//2 + 1
        if self.fitboxsize % 2 == 0:
            # for even fitbox sizes, need to truncate ymax/xmax by 1
            ymax -= 1
            xmax -= 1

        data_stamp = data[ymin:ymax, xmin:xmax]
        self.data_stamp = data_stamp

        # store coordinates of stamp also
        dy_img, dx_img = np.indices(data.shape, dtype=float)
        dy_img -= data_center[1]
        dx_img -= data_center[0]

        dx_data_stamp = dx_img[ymin:ymax, xmin:xmax]
        dy_data_stamp = dy_img[ymin:ymax, xmin:xmax]
        self.data_stamp_RA_offset = -dx_data_stamp
        self.data_stamp_Dec_offset = dy_data_stamp
        self.data_stamp_RA_offset_center = self.data_stamp_RA_offset[0, self.fitboxsize // 2]
        self.data_stamp_Dec_offset_center = self.data_stamp_Dec_offset[self.fitboxsize // 2, 0]

        # generate noise map if necessary
        if noise_map is None:
            # blank map
            noise_stamp = np.zeros(data_stamp.shape)

            # define exclusion around planet.
            distance_from_planet = np.sqrt((dx_img - (xguess - data_center[0]))**2 +
                                           (dy_img - (yguess - data_center[1]))**2)
            # define radial coordinate
            rimg = np.sqrt(dx_img**2 + dy_img**2)

            # calculate noise for each pixel in the data_stamp stamp
            for y_index, x_index in np.ndindex(data_stamp.shape):
                r_pix = np.sqrt(dy_data_stamp[y_index, x_index]**2 + dx_data_stamp[y_index, x_index]**2)
                pixels_for_noise = np.where((np.abs(rimg - r_pix) <= dr/2.) & (distance_from_planet > exclusion_radius))
                noise_stamp[y_index, x_index] = np.nanstd(data[pixels_for_noise])

        else:
            noise_stamp = noise_map[ymin:ymax, xmin:xmax]

        self.noise_map = noise_stamp


    def set_kernel(self, covar, covar_param_guesses, covar_param_labels, include_readnoise=False,
                   read_noise_fraction=0.01):
        """
        Set the Gaussian process kernel used in our astrometric fit

        Args:
            covar: Covariance kernel for GP regression. If string, can be "matern32" or "sqexp"
                    Can also be a function: cov = cov_function(x_indices, y_indices, sigmas, cov_params)
            covar_param_guesses: a list of guesses on the hyperparmeteres (size of N_hyperparams)
            covar_param_labels: a list of strings labelling each covariance parameter
            include_readnoise: if True, part of the noise is a purely diagonal term (i.e. read/photon noise)
            read_noise_fraction: fraction of the total measured noise is read noise (between 0 and 1)

        Returns:

        """
        if isinstance(covar, str):
            if covar.lower() == "matern32":
                self.covar = covars.matern32
            elif covar.lower() == "sqexp":
                self.covar = covars.sq_exp
            else:
                raise ValueError("Covariance matricies currently supported are 'matern32' and 'sqexp'")
        else:
            # this better be a covariance function. We're trusting you
            self.covar = covar

        self.covar_param_guesses = covar_param_guesses
        self.covar_param_labels = covar_param_labels

        if include_readnoise:
            self.include_readnoise = True
            self.covar_param_guesses.append(read_noise_fraction)
            self.covar_param_labels.append(r"K_{\delta}")


    def set_bounds(self, dRA, dDec, df, covar_param_bounds, read_noise_bounds=None):
        """
        Set bounds on Bayesian priors. All paramters can be a 2 element tuple/list/array that specifies
        the lower and upper bounds x_min < x < x_max. Or a single value whose interpretation is specified below
        If you are passing in both lower and upper bounds, both should be in linear scale!
        Args:
            dRA: Distance from initial guess position in pixels. For a single value, this specifies the largest distance
                form the initial guess (i.e. RA_guess - dRA < x < RA_guess + dRA)
            dDec: Same as dRA except with Dec
            df: Flux range. If single value, specifies how many orders of 10 the flux factor can span in one direction
                (i.e. log_10(guess_flux) - df < log_10(guess_flux) < log_10(guess_flux) + df
            covar_param_bounds: Params for covariance matrix. Like df, single value specifies how many orders of
                                magnitude parameter can span. Otherwise, should be a list of 2-elem touples
            read_noise_bounds: Param for read noise term. If single value, specifies how close to 0 it can go
                                based on powers of 10 (i.e. log_10(-read_noise_bound) < read_noise < 1 )

        Returns:

        """
        self.bounds = []

        # x/RA bounds
        if np.size(dRA) == 2:
            self.bounds.append(dRA)
        else:
            self.bounds.append([self.guess_RA_offset - dRA, self.guess_RA_offset + dRA])

        # y/Dec bounds
        if np.size(dDec) == 2:
            self.bounds.append(dDec)
        else:
            self.bounds.append([self.guess_Dec_offset - dDec, self.guess_Dec_offset + dDec])

        # flux bounds
        # need to guess flux if None
        if self.guess_flux is None:
            if self.fm_stamp is not None and self.data_stamp is not None:
                # use the two to scale it and put it in log scale
                self.guess_flux = np.max(self.data_stamp) / np.max(self.fm_stamp)
            else:
                # we haven't read in the data_stamp and FM yet. Assume they're on the same scale
                # should be in log scale
                self.guess_flux = 1
        if np.size(df) == 2:
            self.bounds.append(df)
        else:
            self.bounds.append([self.guess_flux / (10.**df), self.guess_flux * (10**df)])

        # hyperparam bounds
        if np.ndim(covar_param_bounds) == 2:
            for covar_param_bound in covar_param_bounds:
                self.bounds.append(covar_param_bound)
        else:
            # this is a 1-D list, with each param specified by one paramter
            for covar_param_bound, covar_param_guess in zip(covar_param_bounds, self.covar_param_guesses):
                self.bounds.append([covar_param_guess / (10.**covar_param_bound),
                                    covar_param_guess * (10**covar_param_guess)])

        if read_noise_bounds is not None:
        # read noise
            if np.size(read_noise_bounds) == 2:
                self.bounds.append(read_noise_bounds)
            else:
                self.bounds.append([10**read_noise_bounds, 1])


    def fit_astrometry(self, nwalkers=100, nburn=200, nsteps=800, save_chain=True, chain_output="bka-chain.pkl",
                       numthreads=None):
        """
        Run a Bayesian fit of the astrometry using MCMC
        Saves to self.chian

        Args:
            nwalkers: number of walkers
            nburn: numbe of samples of burn-in for each walker
            nsteps: number of samples each walker takes
            save_chain: if True, save the output in a pickled file
            chain_output: filename to output the chain to
            numthreads: number of threads to use

        Returns:

        """
        # create array of initial guesses
        # array of guess RA, Dec, and flux
        # for everything that's not RA/Dec offset, should be converted to log space for MCMC sampling
        init_guess = np.array([self.guess_RA_offset, self.guess_Dec_offset, math.log(self.guess_flux)])
        # append hyperparams for covariance matrix, which also need to be converted to log space
        init_guess = np.append(init_guess, np.log(self.covar_param_guesses))
        # number of dimensions of MCMC fit
        ndim = np.size(init_guess)

        # initialize walkers in a ball around the best fit value
        pos = [init_guess + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

        # prior bounds also need to be put in log space
        sampler_bounds = np.copy(self.bounds)
        sampler_bounds[2:] = np.log(sampler_bounds[2:])

        global lnprob
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(self, sampler_bounds, self.covar),
                                        threads=numthreads)

        # burn in
        print("Running burn in")
        pos, _, _ = sampler.run_mcmc(pos, nburn)
        # reset sampler
        sampler.reset()
         
        # chains should hopefulyl have converged. Now run MCMC
        print("Burn in finished. Now sampling posterior")
        sampler.run_mcmc(pos, nsteps)
        print("MCMC sampler has finished")

        # convert chains in log space back in linear space
        sampler.chain[:,:,2:] = np.exp(sampler.chain[:,:,2:])

        # save state
        self.sampler = sampler

        # save best fit values
        # percentiles has shape [ndims, 3]
        percentiles = np.swapaxes(np.percentile(sampler.flatchain, [16, 50, 84], axis=0), 0, 1)
        self.RA_offset = percentiles[0][1]
        self.RA_offset_1sigma = (percentiles[0][0], percentiles[0][2])
        self.Dec_offset = percentiles[1][1]
        self.Dec_offset_1sigma = (percentiles[1][0], percentiles[1][2])
        self.flux = percentiles[2][1]
        self.flux_1sigma = (percentiles[2][0], percentiles[2][2])
        self.covar_param_bestfits = [thispercentile[1] for thispercentile in percentiles[3:]]
        self.covar_param_1sigma = [(thispercentile[0], thispercentile[2]) for thispercentile in percentiles[3:]]

        if save_chain:
            pickle_file = open(chain_output, 'w')
            pickle.dump(sampler.chain, pickle_file)
            pickle.dump(sampler.lnprobability, pickle_file)
            pickle.dump(sampler.acceptance_fraction, pickle_file)
            #pickle.dump(sampler.acor, pickle_file)
            pickle_file.close()


    def make_corner_plot(self, fig=None):
        """
        Generate a corner plot of the posteriors from the MCMC
        Args:
            fig: if not None, a matplotlib Figure object

        Returns:
            fig: the Figure object. If input fig is None, function will make a new one

        """
        all_labels = [r"x", r"y", r"$\alpha$"]
        all_labels = np.append(all_labels, self.covar_param_labels)

        fig = corner.corner(self.sampler.flatchain, labels=all_labels, quantiles=[0.16, 0.5, 0.84], fig=fig)

        return fig


    def best_fit_and_residuals(self, fig=None):
        """
        Generate a plot of the best fit FM compared with the data_stamp and also the residuals
        Args:
            fig (matplotlib.Figure): if not None, a matplotlib Figure object

        Returns:
            fig (matplotlib.Figure): the Figure object. If input fig is None, function will make a new one

        """
        if fig is None:
            fig = plt.figure(figsize=(12, 4))

        # create best fit FM
        dx = -(self.RA_offset - self.data_stamp_RA_offset_center)
        dy = self.Dec_offset - self.data_stamp_Dec_offset_center

        fm_bestfit = self.flux * sinterp.shift(self.fm_stamp, [dy, dx])
        if self.padding > 0:
            fm_bestfit = fm_bestfit[self.padding:-self.padding, self.padding:-self.padding]

        # make residual map
        residual_map = self.data_stamp - fm_bestfit

        # normalize all images to same scale
        colornorm = matplotlib.colors.Normalize(vmin=np.percentile(self.data_stamp, 0.03),
                                                vmax=np.percentile(self.data_stamp, 99.7))

        # plot the data_stamp
        ax1 = fig.add_subplot(131)
        im1 = ax1.imshow(self.data_stamp, interpolation='nearest', cmap='cubehelix', norm=colornorm)
        ax1.invert_yaxis()
        ax1.set_title("Data")
        ax1.set_xlabel("X (pixels)")
        ax1.set_ylabel("Y (pixels)")

        ax2 = fig.add_subplot(132)
        im2 = ax2.imshow(fm_bestfit, interpolation='nearest', cmap='cubehelix', norm=colornorm)
        ax2.invert_yaxis()
        ax2.set_title("Best-fit Model")
        ax2.set_xlabel("X (pixels)")

        ax3 = fig.add_subplot(133)
        im3 = ax3.imshow(residual_map, interpolation='nearest', cmap='cubehelix', norm=colornorm)
        ax3.invert_yaxis()
        ax3.set_title("Residuals")
        ax3.set_xlabel("X (pixels)")

        fig.subplots_adjust(right=0.82)
        fig.subplots_adjust(hspace=0.4)
        ax_pos = ax3.get_position()

        cbar_ax = fig.add_axes([0.84, ax_pos.y0, 0.02, ax_pos.height])
        cb = fig.colorbar(im1, cax=cbar_ax)
        cb.set_label("Counts (DN)")

        return fig


def lnprior(fitparams, bounds):
    """
    Bayesian prior

    Args:
        fitparams: array of params (size N)

        bounds: array of (N,2) with corresponding lower and upper bound of params
                bounds[i,0] <= fitparams[i] < bounds[i,1]

    Returns:
        prior: 0 if inside bound ranges, -inf if outside

    """
    prior = 0.0

    for param, bound in zip(fitparams, bounds):
        if (param >= bound[1]) | (param < bound[0]):
            prior *= -np.inf

    return prior


def lnlike(fitparams, fma, cov_func):
    """
    Likelihood function
    Args:
        fitparams: array of params (size N). First three are [dRA,dDec,f]. Additional parameters are GP hyperparams
                    dRA,dDec: RA,Dec offsets from star. Also coordianaes in self.data_{RA,Dec}_offset
                    f: flux scale factor to normalizae the flux of the data_stamp to the model
        fma (FMAstrometry): a FMAstrometry object that has been fully set up to run
        cov_func (function): function that given an input [x,y] coordinate array returns the covariance matrix
                  e.g. cov = cov_function(x_indices, y_indices, sigmas, cov_params)

    Returns:
        likeli: log of likelihood function (minus a constant factor)
    """
    dRA_trial, dDec_trial, f_trial, hyperparms_trial = fitparams

    # get trial parameters out of log space
    f_trial = math.exp(f_trial)
    hyperparms_trial = np.exp(hyperparms_trial)

    dx = -(dRA_trial - fma.data_stamp_RA_offset_center)
    dy = dDec_trial - fma.data_stamp_Dec_offset_center

    fm_shifted = sinterp.shift(fma.fm_stamp, [dy, dx])

    if fma.padding > 0:
        fm_shifted = fm_shifted[fma.padding:-fma.padding, fma.padding:-fma.padding]

    diff_ravel = fma.data_stamp.ravel() - f_trial * fm_shifted.ravel()

    cov = cov_func(fma.data_stamp_RA_offset.ravel(), fma.data_stamp_Dec_offset.ravel(), fma.noise_map.ravel(),
                   hyperparms_trial)


    # solve Cov * x = diff for x = Cov^-1 diff. Numerically more stable than inverse
    # to make it faster, we comptue the Cholesky factor and use it to also compute the determinent
    (L_cov, lower_cov) = linalg.cho_factor(cov)
    cov_inv_dot_diff = linalg.cho_solve((L_cov, lower_cov), diff_ravel) # solve Cov x = diff for x
    residuals = diff_ravel.dot(cov_inv_dot_diff)

    # compute log(det(Cov))
    logdet = 2*np.sum(np.log(np.diag(L_cov)))
    constant = logdet

    return -0.5 * (residuals + constant)


def lnprob(fitparams, fma, bounds, cov_func):
    """
    Function to compute the relative posterior probabiltiy. Product of likelihood and prior
    Args:
        fitparams: array of params (size N). First three are [dRA,dDec,f]. Additional parameters are GP hyperparams
                    dRA,dDec: RA,Dec offsets from star. Also coordianaes in self.data_{RA,Dec}_offset
                    f: flux scale factor to normalizae the flux of the data_stamp to the model
        fma: a FMAstrometry object that has been fully set up to run
        bounds: array of (N,2) with corresponding lower and upper bound of params
                bounds[i,0] <= fitparams[i] < bounds[i,1]
        cov_func: function that given an input [x,y] coordinate array returns the covariance matrix
                  e.g. cov = cov_function(x_indices, y_indices, sigmas, cov_params)

    Returns:

    """
    lp = lnprior(fitparams, bounds)

    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(fitparams, fma, cov_func)


