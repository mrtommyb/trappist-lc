from __future__ import division, print_function
import numpy as np
from astropy.stats.funcs import median_absolute_deviation as MAD
from scipy.ndimage import label

from martinsff import martinsff
import extract_lc

import yash_bls
from blssearch import doSearch, plotSearch, get_qf
from photo_test import raw_moment, intertial_axis, plot_bars

from ktransit import FitTransit, LCModel
import ktransit
from planet_params import *

def get_label_im(fluxarr,bg_cut):
    fbg = bg_sub(fluxarr)
    flatimx = np.nanmedian(fbg,axis=0)
    vals = flatimx[np.isfinite(flatimx)].flatten()
    mad_cut = 1.4826 * MAD(vals) * bg_cut

    flatimx[np.isnan(flatimx)] = 0.
    region = np.where(flatimx > mad_cut,1,0)
    lab = label(region)[0]

    #find the central pixel
    imshape = np.shape(flatimx)
    centralpix = [1+imshape[0] // 2,1+imshape[1] // 2]

    regnum = lab[centralpix[0],centralpix[1]]
    labim = np.where(lab == regnum, 1, 0)

    return labim


def bg_sub(fla):
    """
    subtract the background from a series of images
    by assuming the aperture is large enough to be
    predominantly background
    """
    for i in xrange(np.shape(fla)[0]):
        fla[i,:,:] = fla[i,:,:] - np.nanmedian(fla[i,:,:])
    return fla

def optimalAperture(t_time, t_fluxarr, t_quality, qual_cut=False, return_qual=False, toss_resat=False, bg_cut=5, skip=0):
    """
    This routine determines an optimal apertures and outputs the flux (i.e. a light curve) from a TPF.


    Inputs:
    ------------
    t_time = 1D array of 'TIME' from TPF
    t_fluxarr = 1D array of 'FLUX' from TPF
    t_quality = 1D array of 'QUALITY' from TPF
    qual_cut = exclude cadences with a non-zero quality flag; this is False by default
    return_qual = if True then nothing is returned; this is True by default
    toss_resat = exclude cadences where there is a wheel resaturation event; this is True by default
    bg_cut = threshold to find pixels that are bg_cut * MAD above the median
    skip = index of first cadence that should be used in the time series

    Outputs:
    ------------
    time = 1D array of time in BKJD
    lc = 1D array of flux measured in optimal aperture
    xbar = 1D array of x-coordinate of target centroid
    ybar = 1D array of y-coordinate of target centroid
    regnum = integer value of brightest pixel
    lab = 2D array identifying pixels used in aperture

    Usage:
    ------------
    tpf,tpf_hdr = ar.getLongTpf(k2id, campaign, header=True)

    tpf_time = tpf['TIME']
    tpf_flux = tpf['FLUX']
    tpf_quality = tpf['QUALITY']

    time,lc,xbar,ybar,regnum,lab = optimalAperture(tpf_time, tpf_flux, tpf_quality, qual_cut=False, return_qual=False, toss_resat=True, bg_cut=5, skip=0)
    """

    time = t_time[skip:]
    fluxarr = t_fluxarr[skip:]
    quality = t_quality[skip:]

    if qual_cut:
        time = time[quality == 0]
        fluxarr = fluxarr[quality == 0,:,:]
    elif toss_resat:
        # cadences where there is a wheel resaturation event
        time = time[quality != 32800]
        fluxarr = fluxarr[quality != 32800,:,:]

    #remove any nans
    try:
      fluxarr[fluxarr == 0] = np.nan
    except ValueError:
      pass

    #subtract background
    flux_b = bg_sub(fluxarr)

    # create a median image to calculate where the pixels to use are
    flatim = np.nanmedian(flux_b,axis=0)

    #find pixels that are X MAD above the median
    vals = flatim[np.isfinite(flatim)].flatten()
    mad_cut = 1.4826 * MAD(vals) * bg_cut

    flatim[np.isnan(flatim)] = 0.
    region = np.where(flatim > mad_cut,1,0)
    lab = label(region)[0]

    #find the central pixel
    imshape = np.shape(flatim)
    centralpix = [1+imshape[0] // 2,1+imshape[1] // 2]

    #find brightest pix within 9x9 of central pix
    #this assumes target is at center of postage stamp which I think is ok
    centflatim = flatim[centralpix[0]-2:centralpix[0]+2,
                centralpix[1]-2:centralpix[1]+2]
    flatimfix = np.where(np.isfinite(centflatim),centflatim,0)
    brightestpix = np.unravel_index(flatimfix.argmax(), centflatim.shape)
    bpixy, bpixx = brightestpix

    #use all pixels in the postage stamp that are X MAD above the median
    #this identifies location of brightest pixel only
    regnum = lab[centralpix[0]-2+bpixy,centralpix[1]-2+bpixx]

    lc = np.zeros_like(time)
    xbar = np.zeros_like(time)
    ybar = np.zeros_like(time)

    #make a rectangular aperture for the moments thing
    ymin = np.min(np.where(lab == regnum)[0])
    ymax = np.max(np.where(lab == regnum)[0])
    xmin = np.min(np.where(lab == regnum)[1])
    xmax = np.max(np.where(lab == regnum)[1])

    momlims = [ymin,ymax+1,xmin,xmax+1]

    #loop that performs the aperture photometry
    for i,fl in enumerate(flux_b):
        lc[i] = np.sum(fl[lab == regnum])
        #lc[i] = np.sum(fl[np.where(lab == 1)]
        momim = fl[momlims[0]:momlims[1],
                    momlims[2]:momlims[3]]
        momim[~np.isfinite(momim)] == 0.0
        xbar[i], ybar[i], cov = intertial_axis(momim)


    if return_qual:
        return None
    else:
        # TODO: think about whether this should be normalized
        return (time,lc, xbar - np.nanmean(xbar), ybar - np.nanmean(ybar), regnum, lab)


def get_lc(time1, fluxarr1, quality1, n_chunks, bg_cut, flatlc_window, smooth_window):


    time, lc, xbar, ybar, regnum, lab = optimalAperture(
            time1, fluxarr1, quality1, qual_cut=False, return_qual=False,
                       toss_resat=False, bg_cut=bg_cut,  skip=None)

    m1 = np.isfinite(lc) * np.isfinite(lc)

    time = time1[m1]
    lc = lc[m1]
    xbar = xbar[m1]
    ybar = ybar[m1]



    flatlc = extract_lc.medfilt(time,lc,window=flatlc_window)

    cadstep = np.int(np.floor(len(time) / n_chunks)) #600
    zpt = len(time) % cadstep
    if zpt==cadstep:
        zpt = 0

    outflux, correction, thr_cad = extract_lc.run_C0_detrend(
            time, flatlc, xbar, ybar, cadstep=cadstep, skip=None)

    not_thr = ~thr_cad
    corflux = (lc[zpt:][not_thr]/
        np.median(lc[zpt:][not_thr])/
        correction[not_thr])

    corflatflux = (flatlc[zpt:][not_thr]/
        np.median(flatlc[zpt:][not_thr])/
        correction[not_thr])

    # The 1.4826 and *4 factors make this similar to a 4-sigma cut.
    mad_cut = 1.4826*MAD(corflatflux-1.) * 4.0
    keep = np.abs(corflatflux-1.) < mad_cut # this might not be used

    m2 = np.ones_like(corflatflux,dtype=bool)#corflatflux < 1.1
    t1 = time[zpt:][not_thr][m2]
    cfflux = extract_lc.medfilt(t1,corflux,window=smooth_window) - 1.0

    return time, lc, xbar, ybar, t1, corflux, cfflux


def transit_fit(t1, cfflux):
    addtime = 4833 # convert from trappist time to kepler time

    time = t1 [(cfflux < 0.005) * (cfflux > -0.015)]  +addtime        # you need a time and a flux
    flux = cfflux[(cfflux < 0.005) * (cfflux > -0.015)]             # there are no transits here :(
    ferr = np.ones_like(time) * 0.001     # uncertainty on the data

    fitT = FitTransit()
    fitT.add_guess_star(rho=50.0, ld1 = 0.43, ld2 = 0.14)  # fixed because I'm sleepy

    for planet in [planetb, planetc, planetd, planete, planetf, planetg, planeth]:
        fitT.add_guess_planet(
            period=planet['period_days'][0], impact=planet['impact'][0],
            T0=planet['t0'][0], rprs=(planet['td_percent'][0]/100)**0.5)

    fitT.add_data(time=time, flux=flux, ferr=ferr, itime=np.ones_like(time) * 0.0188)

    vary_star = ['rho', 'zpt']      # free stellar parameters
    vary_planet = (['period',       # free planetary parameters
            'T0', 'impact',
            'rprs'])                # free planet parameters are the same for every planet you model

    fitT.free_parameters(vary_star, vary_planet)
    fitT.do_fit()                   # run the fitting

    return time, flux, fitT

def get_qf(time,flux,epoch,period,transitmodel=None):
    date1 = (time - epoch) + 0.5*period
    phi1 = (((date1 / period) - np.floor(date1/period)) * period) - 0.5*period
    q1 = np.sort(phi1)
    f1 = (flux[np.argsort(phi1)]) * -1.E6
    if transitmodel is not None:
        m1 = (transitmodel[np.argsort(phi1)]) * -1.E6
        return q1,f1,m1
    else:
        return q1,f1

def bin_data(phi,flux,bins,model=None):
    phi = np.array(phi)
    flux = np.array(flux)
    phibin = []
    fluxbin = []
    stdbin = []
    bins=int(bins)
    for i in (bins*np.arange(len(phi)//bins))+(bins//2):
        if model == None:
            goodpoints = np.ones(len(flux),dtype=bool)
        else:
            goodpoints = flux-model < 3* np.std(flux-model)
        flux2 = flux[goodpoints]
        phi2 = phi[goodpoints]
        phibin.append(np.median(phi2[i-bins//2:i+bins//2]))
        fluxbin.append(np.median(flux2[i-bins//2:i+bins//2]))
        stdbin.append(np.std(flux2[i-bins//2:i+bins//2]))
    return np.array(phibin), np.array(fluxbin), np.array(stdbin) / np.sqrt(bins)


def model_ktransit(time, rho, zpt, ld1, ld2, fitresultsplanet):
    M = LCModel()
    M.add_star(rho=rho,zpt=zpt,ld1=ld1,
            ld2=ld2,dil=0)
    M.add_planet(T0=fitresultsplanet['T0'],
                 period=fitresultsplanet['period'],
                 impact=fitresultsplanet['impact'],
                 rprs=fitresultsplanet['rprs']
                 )
    M.add_data(time=time)
    return M.transitmodel

def plot_model(ax,time,flux, rho, zpt, ld1, ld2, fitresultsplanet):

    T0=fitresultsplanet['T0']
    period=fitresultsplanet['period']
    q1, f1 = get_qf(time,flux,T0,period)

    model = model_ktransit(time, rho, zpt, ld1, ld2, fitresultsplanet)
    q1, m1 = get_qf(time,model,T0,period)
    ax.scatter(q1*24., f1/10000,color='k',alpha=period**0.5/5,s=4
              )
    ax.plot(q1*24., m1/10000,color='r')
    bq, bf, be = bin_data(q1,f1,bins=50.//period)
    ax.errorbar(bq*24., bf/10000, yerr=be/10000, ls='',color='b')
    ax.set_xlim([-4,4])
    ax.set_ylim([1.19,-0.75])
    #ax.set_xlabel('Time from mid-transit (days)')
    ax.set_ylabel('Transit depth (%)', fontsize=17)
    ax.minorticks_on()
    return ax

