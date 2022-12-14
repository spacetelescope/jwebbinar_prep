# !/usr/bin/env python3.6
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
from scipy.optimize import nnls
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.table import Table
import time

# Fundamental constants
k_b = 1.38065e-23  # m^2 kg s^-2 K^-1 (SI)
c = 299792458.  # m*sec^-1
c_microns = 299792458000000  # microns*sec^-1
h = 6.62607e-34  # m^2 kg s^-1


def blackbody(wave, temperature):
    """
    Black body function, units are in SI
    :param wave: array of floats
        wavelength array, steps must be regular
    :param temperature: float
        black body temperature (in K)
    :return black body luminosity profile
    """
    luminosity = (2 * h * (c ** 2)) / ((wave ** 5) * (np.exp(h * c / (wave * k_b * temperature)) - 1))
    return luminosity


def gauss(wave, peak_center, sigma):
    """
    returns Gaussian function on the wavelength array given
    :param wave: array
        wavelengths array
    :param peak_center: float
        position of the center of the peak
    :param sigma: float
        width of the gaussian
    :return array of gaussian values centered at the given peak
    """
    return np.exp(-np.power(wave - peak_center, 2.) / (2 * np.power(sigma, 2.)))


class Spectrum:
    """
    Use a model to fit observed spectra to apply a new NMF on real data
    Needs to initialize with some parameters first:
        the spectrum to fit, the corresponding wavelength array and the spectral resolution of the spectrum to fit
    """

    def __init__(self, spectrum_to_fit, wave_axis, resolution, redshift=0., templates='Foschino', first_run=True):
        """
        Initialize the class using the following parameters:
        :param spectrum_to_fit: array
            The spectrum to fit
        :param wave_axis: array
            The array of wavelengths corresponding to the spectrum to fit, in microns
        :param resolution: float
            The spectrum to fit resolution corresponding to R=lambda/delta_lambda
        :param redshift: float
            The redshift of the source spectrum, default is 0.
        :param templates: string
            Name of the PAHTAH templates to use, default is 'Foschino'
        :param first_run: bool
            default is True, if false, means it will load the files gas_lines.fits, continuum.fits and AIBs.fits
            to have a quicker fit

        """
        if templates == 'Foschino':
            wave_elementary_axis_filename = '/home/shared/preloaded-fits/pahtat/Elementary_spectra_Foschino_wave.fits'
            elementary_spectra = '/home/shared/preloaded-fits/pahtat/Elementary_spectra_Foschino.fits'
        else:
            raise NameError('Templates name is not a valid one, try to type "Foschino"')
        self.spectrum_to_fit = spectrum_to_fit
        self.redshift = redshift
        self.wave = wave_axis / (1 + self.redshift)
        self.wave_elementary = fits.getdata(wave_elementary_axis_filename)
        self.wave_elementary *= 1e6  # because we need wavelengths in microns
        self.elementary_spectra = fits.getdata(elementary_spectra)
        self.wave_resolution = resolution
        first_x_value = max(self.wave[0], self.wave_elementary[0])
        last_x_value = min(self.wave[-1], self.wave_elementary[-1])
        delta_lambda = first_x_value / self.wave_resolution
        x_to_fit = np.arange(first_x_value, last_x_value, delta_lambda)

        if self.wave_elementary.shape != self.wave.shape:
            interpolate_elementary = interpolate.interp1d(self.wave_elementary, self.elementary_spectra)
            interpolate_to_fit = interpolate.interp1d(self.wave, self.spectrum_to_fit)
            self.spectrum_to_fit = interpolate_to_fit(x_to_fit)
            self.elementary_spectra = interpolate_elementary(x_to_fit)
            self.wave = x_to_fit

        if first_run is False:
            self.continuum = fits.getdata('Data/continuum.fits')
            self.gas_lines = fits.getdata('Data/gas_lines.fits')
            self.AIBs = fits.getdata('Data/AIBs.fits')
        else:
            self.continuum = None
            self.gas_lines = None
            self.AIBs = np.zeros([self.wave.shape[0], self.elementary_spectra.shape[0]])
            self.PAHp = np.zeros([self.wave.shape[0], 1])
            self.PAH0 = np.zeros([self.wave.shape[0], 1])
            self.eVSG = np.zeros([self.wave.shape[0], 1])
            self.PAHx = np.zeros([self.wave.shape[0], 1])

        self.AIBs_gauss = None
        self.model = None
        self.fit_result = None
        self.nb_comp = None
        self.c_ext = None
        self.gas_lines_ind = None
        self.species = None
        self.central_lambda = None
        self.first_run = first_run

    def build_continuum(self, t_min=40, t_max=6000, delta_t=20, nh=None):
        """
        Creates the continuum model
        :param t_min: float
            Minimal temperature of the blackbodies catalog in K, default is 40
        :param t_max: float
            Maximal temperature of the blackbodies catalog in K, default is 6000
        :param delta_t: float
            step size in temperature between two blackbodies, in K, default is 20
        :param nh: array
            column densities taken into account for the extinction, default is None

        fills the continuum Class parameter with an array of nb_CN columns of len(wave) points
        Writes the continuum in 'Data/continuum.fits' file
        """

        # extinction law from Weingartner and Drain 2001, Rv=5.5
        if nh is None:
            nh = [0, 1e18, 1.0e19, 1.0e20, 1.6e22, 5e22, 1e23]
        wave_draine__, albe, cos, cext__ = np.loadtxt('Data/extinction_curve/Cext_tabulation_draine_55.txt',
                                                      unpack=True,
                                                      usecols=[0, 1, 2, 3])
        Cext_ = cext__[np.where(wave_draine__ < (self.wave[-1]))]
        wave_draine_ = wave_draine__[np.where(wave_draine__ < (self.wave[-1]))]

        cext_ = Cext_[np.where(wave_draine_ > (self.wave[0]))]
        wave_draine = wave_draine_[np.where(wave_draine_ > (self.wave[0]))]

        C_ext_new = interpolate.interp1d(wave_draine, cext_, fill_value="extrapolate")
        Cext = C_ext_new(self.wave)
        self.c_ext = Cext
        # =============================================================================
        # Grey body parameters
        # number of black body
        nb_CN = len(nh) * ((t_max - t_min) / delta_t)
        nb_CN = int(nb_CN)
        slope = np.arange(0.1, 3, 0.1)
        # creates continuum model
        self.continuum = np.zeros([len(self.wave), nb_CN + slope.shape[0]])
        p = 0
        for num_discr_cn in range(len(nh)):
            extinction = np.exp(-Cext * nh[num_discr_cn])
            for cn in range(int(nb_CN / len(nh))):
                # grey body
                grey_body = extinction * blackbody(self.wave * 1e-6, t_min + cn * delta_t)
                integral = np.sum(grey_body)
                self.continuum[:, p] = grey_body / integral
                p += 1
        for factor in slope:
            self.continuum[:, p] = factor * self.wave
            p += 1

        print('continuum array shape:', self.continuum.shape)
        fits.writeto('Data/continuum.fits', self.continuum, overwrite=True)

    def build_gaslines(self):
        """
        Function generating the gas lines model.

        Fills the gas lines fitting catalog into the gas-lines parameter, an array of catalog Gaussian columns and
        len(wave) spectral points
        Write the gas lines in 'Data/gas_lines.fits' file
        """
        table_lines = Table.read('Data/line_list_all.txt', format='ascii')
        central_lambda = table_lines['wavelength']
        self.species = table_lines['Species']
        self.central_lambda = central_lambda

        width = (central_lambda / self.wave_resolution)
        self.gas_lines = np.zeros([len(self.wave), np.size(central_lambda)])

        for i in range(len(central_lambda)):
            self.gas_lines[:, i] = gauss(self.wave, (central_lambda[i]), width[i] / 2.355)

        print('gas array shape:', self.gas_lines.shape)
        fits.writeto('Data/gas_lines.fits', self.gas_lines, overwrite=True)

    def build_aib(self, norm):
        """
        Function generating the catalog used to fit the AIBs (Aromatic Infrared Band).
        :param norm: bool
            choose or not to normalize the elementary spectra if needed.
        Fills the AIB parameter with fitting catalog, an array of 4 columns and len(wave) spectral points
        Writes the AIB fitting catalog into 'Data/AIBs.fits' file
        """
        if norm:
            for i_norm in range(self.elementary_spectra.shape[0]):
                print(self.elementary_spectra.shape, self.wave.shape)
                self.AIBs[:, i_norm] = self.elementary_spectra[i_norm] / np.trapz(self.elementary_spectra[i_norm],
                                                                                  self.wave)
            self.PAHp[:, 0] = self.AIBs[:, 0]
            self.PAH0[:, 0] = self.AIBs[:, 2]
            self.eVSG[:, 0] = self.AIBs[:, 1]
            self.PAHx[:, 0] = self.AIBs[:, 3]
        else:
            for i_norm in range(self.elementary_spectra.shape[0]):
                print(self.elementary_spectra.shape, self.wave.shape)
                self.AIBs[:, i_norm] = self.elementary_spectra[i_norm].copy

        print('AIBs arrays shape:', self.AIBs.shape)
        fits.writeto('Data/AIBs.fits', self.AIBs, overwrite=True)

    def build_catalog(self):
        """
        Construction of the basis of the model
        """
        print('Fit of the observed spectra with the elementary spectra for the AIBs')
        if self.first_run:
            self.build_aib(norm=True)
            self.build_gaslines()
            self.build_continuum()

        nb_comp = [self.gas_lines.shape[1], self.continuum.shape[1], self.PAHp.shape[1], self.PAH0.shape[1],
                   self.eVSG.shape[1], self.PAHx.shape[1]]
        self.model = np.concatenate((self.gas_lines, self.continuum, self.PAHp, self.PAH0, self.eVSG, self.PAHx), 1)
        self.nb_comp = nb_comp

    def fit(self, fit=True):
        """
        Fitting function using a non-negative least square algorithm, calculating the scale factor of each component
        of the base.
        It uses the "build_catalog" method the first time if no catalog has been generated.
        fill the "fit_result" object with a dictionary containing the fitting information
        :param fit: bool
            Runs PAHTAT on the spectrum to fit. Default is True.
        """
        if fit:
            self.build_catalog()
            self.fit_result = {}

        obj_total_flux = np.zeros([self.spectrum_to_fit.shape[0], 1])
        # fitting
        print('Fitting spectrum...')
        start_time = time.time()
        x, norm1 = nnls(self.model, self.spectrum_to_fit)
        print("--- %s seconds ---" % (time.time() - start_time))
        print(norm1)
        # extracting each component
        model = np.sum(x * self.model, axis=1)
        save_x = x.copy()
        obj_total_flux[0] = np.trapz(self.spectrum_to_fit, self.wave)

        if len(self.nb_comp) > 1:
            gas = np.sum(x[:self.nb_comp[0]] * self.model[:, :self.nb_comp[0]], axis=1)
            self.gas_lines_ind = x[:self.nb_comp[0]] * self.model[:, :self.nb_comp[0]]

            cont = np.sum(x[self.nb_comp[0]:(self.nb_comp[0] + self.nb_comp[1])] *
                          self.model[:, self.nb_comp[0]:(self.nb_comp[0] + self.nb_comp[1])], axis=1)

            AIBs = np.sum(x[(self.nb_comp[0] + self.nb_comp[1]):] *
                          self.model[:, (self.nb_comp[0] + self.nb_comp[1]):], axis=1)

            pahp_position_l = self.nb_comp[0] + self.nb_comp[1]
            pahp_position_h = self.nb_comp[0] + self.nb_comp[1] + self.nb_comp[2]
            pah0_position = self.nb_comp[0] + self.nb_comp[1] + self.nb_comp[2] + self.nb_comp[3]
            evsg_position = self.nb_comp[0] + self.nb_comp[1] + self.nb_comp[2] + self.nb_comp[3] + self.nb_comp[4]
            pahx_position = self.nb_comp[0] + self.nb_comp[1] + self.nb_comp[2] + self.nb_comp[3] + self.nb_comp[4] \
                            + self.nb_comp[5]
            PAHp = np.sum(x[pahp_position_l:pahp_position_h] * self.model[:, pahp_position_l:pahp_position_h], axis=1)
            PAH0 = np.sum(x[pahp_position_h:pah0_position] * self.model[:, pahp_position_h:pah0_position], axis=1)
            eVSG = np.sum(x[pah0_position:evsg_position] * self.model[:, pah0_position:evsg_position], axis=1)
            PAHx = np.sum(x[evsg_position:pahx_position] * self.model[:, evsg_position:pahx_position], axis=1)

            self.fit_result = {'wave': self.wave, 'model': model, 'AIBs': AIBs, 'continuum': cont, 'gas': gas,
                               'Base': self.model, 'weights': save_x, 'int': obj_total_flux, 'PAHp': PAHp, 'PAH0': PAH0,
                               'eVSG': eVSG, 'PAHx': PAHx}

        else:
            self.fit_result = {'wave': self.wave, 'model': model, 'Base': self.model, 'weights': save_x,
                               'int': obj_total_flux}

    def plot_results(self, fit=True):
        """
        Plots the fit of the spectrum
        :param fit: bool
            if True, fits the spectrum using pahtat templates, default is True

        """
        if fit:
            self.fit()

        plt.figure(3, figsize=(16, 6))
        plt.plot(self.wave, self.spectrum_to_fit, 'blue', label='data')
        plt.plot(self.wave, self.fit_result['continuum'], 'black', ls='--', label='continuum')
        plt.plot(self.wave, self.fit_result['gas'], 'k', label='gas')
        plt.plot(self.wave, self.fit_result['model'], 'red', label='fit', zorder=9)
        plt.plot(self.wave, self.fit_result['PAHp'], label=r'PAH$^+$')
        plt.plot(self.wave, self.fit_result['PAH0'], label=r'PAH$^0$')
        plt.plot(self.wave, self.fit_result['eVSG'], label=r'eVSG')
        plt.plot(self.wave, self.fit_result['PAHx'], label=r'PAH$^x$')

        plt.legend(fontsize=15)
        plt.xlabel(r"Wavelength ($\mu$m)", fontsize=20)
        plt.ylabel(r'Flux density (MJy/sr)', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        # plt.subplots_adjust(top=0.989, bottom=0.175, left=0.08, right=0.995, hspace=0.2, wspace=0.2)
        plt.tight_layout()
        plt.show()

    def save_results(self, filename, gas_lines=False):
        """
        saves the fitting results in an ascii file with filename.txt
        """
        results = Table()
        results['wavelength'] = self.wave
        results['data'] = self.spectrum_to_fit
        results['gas_lines'] = self.fit_result['gas']
        results['continuum'] = self.fit_result['continuum']

        results['AIB_1'] = self.fit_result['PAHp']
        results['AIB_2'] = self.fit_result['eVSG']
        results['AIB_3'] = self.fit_result['PAH0']
        results['AIB_4'] = self.fit_result['PAHx']
        results['model'] = self.fit_result['model']
        results.write(filename + '.txt', format='ascii', overwrite=True)

        if gas_lines:
            gas_contrib = Table()
            gas_contrib['species'] = self.species
            gas_contrib['central_lambda'] = self.central_lambda
            flux = []
            for i in range(np.size(self.species)):
                # This computes the line fluxes in erg/s/cm2 assuming the spectrum is in MJy/sr
                flux.append(np.abs(np.trapz(np.transpose(self.gas_lines_ind)[i] * 1e-17, c / (self.wave * 1e-6))))
            gas_contrib['flux'] = flux
            gas_contrib.write(filename + '_line_fluxes.txt', format='ascii', overwrite=True)
