# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

from os.path import join
import numpy as np
from .._in import get_in_params as g
from get_met_ages_values import get_m_a_vls as gmav
from get_isochs import get_isochs as gi
from cmd_phot_systs import phot_wavelengths as pw
from get_CCM_coefs import ccm_model as gcc


def get_mag_idx(mi, phot_params):
    '''
    Return index of magnitude in stored isochrones columns.
    '''
    # Start on 2 to account for the initial and actual mass columns.
    main_mag_idx = 2
    for sys in phot_params[2]:
        for mag in sys[1]:
            if mag == mi:
                # Index of magnitude in stored theoretical values obtained.
                # Get wavelength in inverse microns for this magnitude.
                mw = pw(sys[0], mag)
                # Get CCM coeficient for this magnitude.
                ccm_c = gcc(mw)
                return main_mag_idx, ccm_c
            else:
                main_mag_idx += 1


def get_mc_order(phot_params, isochs_interp):
    '''
    Order the magnitudes and colors in the same way as those stored in the
    input photometric data file.
    Obtain the CCM extinction coeficients for each magnitude and color.
    '''

    isoch_order, ccm_coefs = [], [[], []]
    for _met in isochs_interp:
        met = []
        for _isoch in _met:
            # Append masses and make room for mags and colors.
            isoch = [_isoch[0], _isoch[1], [], []]
            # Search and append ordered magnitudes.
            for m in phot_params[1][0]:
                # Return index of this mag as stored in phot_params[2]
                mi, ccm_i = get_mag_idx(m[1:], phot_params)
                isoch[2].append(_isoch[mi])
                ccm_coefs[0].append(ccm_i)
            # Search and append ordered colors.
            for c in phot_params[1][1]:
                # Return index of this mag as stored in phot_params[2]
                c1, c2 = c[1:].split('-')
                # Search for indexes of each magnitude in the color.
                m1, ccm_1 = get_mag_idx(c1, phot_params)
                m2, ccm_2 = get_mag_idx(c2, phot_params)
                # Append color.
                isoch[3].append((_isoch[m1] - _isoch[m2]))
                ccm_coefs[1].append((ccm_1 - ccm_2))
            met.append(isoch)
        isoch_order.append(met)

    return isoch_order, ccm_coefs


def interp_isoch(isoch_list):
    '''
    Interpolate extra points into the theoretical isochrones passed.
    '''
    N = 2000
    t = np.linspace(0, 1, N)

    isochs_interp = [[] for _ in isoch_list]
    # For each metallicity value.
    for i, met in enumerate(isoch_list):
        # For each age value (ie: isochrone)
        for isoch in met:
            # Interpolate all magnitudes defined.
            xp = np.linspace(0, 1, len(isoch[0]))
            # Store isochrone's interpolated values.
            isochs_interp[i].append(np.asarray([np.interp(t, xp, _)
                for _ in isoch]))

    # *** DELETE ****
    #isoch_inter = np.asarray([np.array(_) for _ in isochrone])

    return isochs_interp


def ip(mypath, phot_params):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''

    ip_list = []
    # Only read files of best fit method is set to run.
    bf_flag = g.bf_params[0]
    if bf_flag is True:

        # Obtain allowed metallicities and ages. Use the first photometric
        # system defined.
        # *WE ASUME ALL PHOTOMETRIC SYSTEMS CONTAIN THE SAME NUMBER OF
        # METALLICITY FILES*
        iso_select = g.ps_params[0]
        iso_path = join(mypath + '/isochrones/' + iso_select + '_' +
            phot_params[2][0][0])
        param_ranges, param_rs, met_f_filter, met_values, age_values = \
        gmav(iso_path)

        print '\n', 'met_vals', met_values
        print 'age_vals', age_values, '\n'

        print 'Interpolating all isochrones for each photometric system.\n'
        # Get isochrones for every photometric system defined.
        isochs_interp = []
        for sys_idx, syst in enumerate(phot_params[2]):

        # *WE ASUME ALL ISOCHRONES OF EQUAL METALLICITY HAVE THE SAME NUMBER
        # OF MASS VALUES ACROSS PHOTOMETRIC SYSTEMS*

            # Get isochrones and their parameter values.
            isoch_list = gi(mypath, met_f_filter, age_values, syst, sys_idx)
            # isoch_list = [met_1, met_2, ..., met_P]

            # Interpolate extra points into all isochrones.
            isochs_interp0 = interp_isoch(isoch_list)
            # For the first photometric system, store the masses that were
            # read.
            if sys_idx == 0:
                isochs_interp.extend(isochs_interp0)
            # For the rest of the systems, just append its magnitudes for
            # each isochrone since we assume the masses are equal.
            else:
                for m_i, _m in enumerate(isochs_interp0):
                    for a_i, _a in enumerate(_m):
                        for mag in _a:
                            isochs_interp[m_i][a_i] = np.append(
                                isochs_interp[m_i][a_i], [mag], 0)

        # isochs_interp = [metal_1, ..., metal_P]
        # metal_i =[age_i, ..., age_Q]
        # age_i = [mass_i, mass_a, mag1, ..., mag_N]
        #
        # mag_1, ..., mag_N are the magnitudes defined in *all* the photom
        # systems, stored in the order presented in photom_params[2]

        # Generate colors (if any) and order magnitudes/colors to match the
        # order of the input photometric data.
        isochs_order, ccm_coefs = get_mc_order(phot_params, isochs_interp)
        # isochs_order = [metal_1, ..., metal_P]
        # metal_i =[age_i, ..., age_Q]
        # age_i = [mass_i, mass_a, [mag1, ..., magN], [col1, ..., colM]

        # Pack params.
        param_values = [met_values, age_values] + param_ranges[2:]
        ip_list = [isochs_order, ccm_coefs, param_values, param_rs]

        iso_ver = {'10': '1.0', '11': '1.1', '12': '1.2S'}
        print ("PARSEC v{} theoretical isochrones read,".format(
            iso_ver[iso_select[-2:]]))
        lens = [len(_) for _ in param_values]
        total = reduce(lambda x, y: x * y, lens, 1)
        print ("interpolated and stored:\n"
        "  {} metallicity values (z),\n"
        "  {} age values (per z),\n"
        "  {} reddening values,\n"
        "  {} distance values,\n"
        "  {} mass values,\n"
        "  {} binary fraction values.".format(*lens))
        print "  = {:.1e} approx total models.".format(total)

    return ip_list