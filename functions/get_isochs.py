# -*- coding: utf-8 -*-
"""
Created on Tue Dic 2 12:00:00 2014

@author: gabriel
"""

import re
import numpy as np
import girardi_isochs_format as gif


def read_met_file(met_f, age_values, line_start, mass_i, mass_a, mags_idx,
    age_format, sys_idx):
    '''
    Read a given metallicity file and return the isochrones for the ages
    within the age range.
    '''

    # Initialize list that will hold all the isochrones for this
    # metallicity value.
    met_i = []

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:

        # Initial empty list for masses and magnitudes.
        isoch_mas_mag = []
        # Initial value for age to avoid 'not defined' error.
        age = -99.

        # Iterate through each line in the file.
        for line in f_iso:

            # Identify beginning of a defined isochrone.
            if line.startswith(line_start):

                # Save stored values if these exist.
                # Skip first age for which the lists will be empty.
                if isoch_mas_mag:
                    # Store magnitudes and masses for this isochrone.
                    met_i.append(isoch_mas_mag)
                    # Reset list.
                    isoch_mas_mag = []

                # Read age value for this isochrone.
                age0 = re.findall(age_format, line)  # Find age in line.
                age = np.around(np.log10(float(age0[0])), 2)

            # If age value falls inside the given range, store the
            # isochrone's data.
            if age in age_values:

                # Save mag, color and mass values for each isochrone.
                if not line.startswith("#"):
                    # Split line.
                    reader = line.split()
                    # Only save masses for one photometric system since
                    # its values are equivalent for the same metallicty and
                    # age across photometric systems.
                    if sys_idx == 0:
                        # Store masses.
                        isoch_mas_mag.append(float(reader[mass_i]))
                        isoch_mas_mag.append(float(reader[mass_a]))
                    # Store defined magnitudes.
                    for mag_i in mags_idx:
                        isoch_mas_mag.append(float(reader[mag_i]))

        # Save the last isochrone when EOF is reached.
        else:
            # If list is not empty.
            if isoch_mas_mag:
                # Store masses and magnitudes for this isochrone.
                met_i.append(isoch_mas_mag)

    # Number of columns stored. Add 2 to account for the 2 masses *only*
    # if this is the first run, ie: for the first photom system. The rest
    # of the phot systems (if more are used) do not have their masses read
    # since mass values are equivalent for equal metallicities and ages.
    n_cols = 2 + len(mags_idx) if sys_idx == 0 else len(mags_idx)
    # Split read values so that the resulting list looks like this:
    #
    # metal_isochs = [age_1, age_2, ..., age_N]
    # age_i = [mass_i, mass_a, mag_1, mag_2, ..., _mag_M]
    #
    # Each mass_x and mag_x are lists that hold all theoretical values
    # read from the metallicity file for that given age.
    metal_isochs = []
    for a_i in range(len(age_values)):
        met_i0 = np.split(np.asarray(met_i[a_i]), len(met_i[a_i]) / n_cols)
        metal_isochs.append(zip(*met_i0))

    return metal_isochs


def get_isochs(mypath, met_f_filter, age_values, syst, sys_idx):
    '''
    Stores the available isochrones of different metallicities and
    ages, according to the ranges given to these parameters.
    '''

    # Read line start format and columns indexes for the selected set of
    # Girardi isochrones.
    line_start, mass_i, mass_a, mags_idx = gif.i_format(syst)
    age_format = gif.age_f()

    # Lists that store the masses and magnitudes of each isochrone.
    # isoch_list = [metal_1, ..., metal_M]
    # metal_i = [isoch_i1, ..., isoch_iN]
    # isoch_ij = [mass_i, mass_a, mag1, mag2, ..., magM]
    # isoch_list[i][j] --> i: metallicity index ; j: age index
    isoch_list = []

    # Iterate in order through all the metallicity files stored for the
    # selected set of isochrones.
    for met_file in met_f_filter:

        metal_isoch = read_met_file(met_file, age_values, line_start, mass_i,
            mass_a, mags_idx, age_format, sys_idx)

        # Store list holding all the isochrones with the same metallicity
        # in the final isochrone list.
        isoch_list.append(metal_isoch)

    return isoch_list