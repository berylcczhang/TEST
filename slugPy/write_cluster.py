"""
Function to write a set of output cluster files in SLUG2 format,
starting from a cluster data set as returned by read_cluster. This can
be used to translate data from one format to another (e.g., bin to
fits), or to consolidate multiple runs into a single output file.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")


def write_cluster(data, model_name, fmt):
    """
    Function to write a set of output cluster files in SLUG2 format,
    starting from a cluster data set as returned by read_cluster.

    Parameters
    ----------
    data : namedtuple
       Cluster data to be written, in the namedtuple format returned
       by read_cluster
    model_name : string
       Base file name to give the model to be written. Can include a
       directory specification if desired.
    fmt : string
       Format for the output file. Allowed values are 'ascii', 'bin'
       or 'binary, and 'fits'.

    Returns
    -------
    Nothing
    """

    # Make sure fmt is valid
    if fmt != 'ascii' and fmt != 'bin' and fmt != 'binary' and \
       fmt != 'fits':
        raise ValueError("fmt must be ascii, bin, binary, or fits")

    # Make sure we're not trying to do fits if we don't have astropy
    if fmt == 'fits' and fits is None:
        raise ValueError("Couldn't import astropy, so fits format "+
                         "is unavailable.")

    ################################################################
    # Write the properties file if we have the data for it
    ################################################################
    if 'form_time' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_prop.txt', 'w')

            # Write header lines
            fp.write(("{:<14s}"*9).
                     format('UniqueID', 'Time', 'FormTime',
                            'Lifetime', 'TargetMass',
                            'BirthMass', 'LiveMass',
                            'NumStar', 'MaxStarMass') + "\n")
            fp.write(("{:<14s}"*9).
                     format('', '(yr)', '(yr)',
                            '(yr)', '(Msun)',
                            '(Msun)', '(Msun)',
                            '', '(Msun)') + "\n")
            fp.write(("{:<14s}"*9).
                     format('-----------', '-----------', '-----------',
                            '-----------', '-----------',
                            '-----------', '-----------',
                            '-----------', '-----------') + "\n")

            # Write data
            for i in range(len(data.id)):
                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*(9*14-3)+"\n")
                fp.write("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   "
                         "{:11.5e}   {:11.5e}   {:11.5e}   {:11d}   "
                         "{:11.5e}\n".
                         format(data.id[i], data.time[i], 
                                data.form_time[i], data.lifetime[i],
                                data.target_mass[i],
                                data.actual_mass[i],
                                data.live_mass[i],
                                data.num_star[i],
                                data.max_star_mass[i]))

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_prop.bin', 'wb')

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_end = ptr + np.argmin(
                        np.logical_and(data.trial[ptr:] == trials[i],
                                       data.time[ptr:] == times[j]))

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.form_time[k])
                        fp.write(data.lifetime[k])
                        fp.write(data.target_mass[k])
                        fp.write(data.actual_mass[k])
                        fp.write(data.live_mass[k])
                        fp.write(data.num_star[k])
                        fp.write(data.max_star_mass[k])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
            cols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
            cols.append(fits.Column(name="FormTime", format="1D",
                                    unit="yr", array=data.form_time))
            cols.append(fits.Column(name="Lifetime", format="1D",
                                    unit="yr", array=data.lifetime))
            cols.append(fits.Column(name="TargetMass", format="1D",
                                    unit="Msun", array=data.target_mass))
            cols.append(fits.Column(name="BirthMass", format="1D",
                                    unit="Msun", array=data.actual_mass))
            cols.append(fits.Column(name="LiveMass", format="1D",
                                    unit="Msun", array=data.live_mass))
            cols.append(fits.Column(name="NumStar", format="1K",
                                    unit="", array=data.num_star))
            cols.append(fits.Column(name="MaxStarMass", format="1D",
                                    unit="Msun", array=data.max_star_mass))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_prop.fits',
                            clobber=True)


    ################################################################
    # Write spectra file if we have the data for it
    ################################################################
    if 'spec' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_spec.txt', 'w')

            # Write header lines
            fp.write(("{:<14s}"*4).
                     format('UniqueID', 'Time', 'Wavelength',
                            'L_lambda') + "\n")
            fp.write(("{:<14s}"*4).
                     format('', '(yr)', '(Angstrom)', '(erg/s/A)')
                     + "\n")
            fp.write(("{:<14s}"*4).
                     format('-----------', '-----------', '-----------',
                            '-----------') + "\n")

            # Write data
            for i in range(len(data.id)):
                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*(9*14-3)+"\n")
                for j in range(len(data.wl)):
                    fp.write("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}"
                             .format(data.id[i], data.time[i],
                                     data.wl[j], data.spec[i,j]) 
                             + "\n")
            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_spec.bin', 'wb')

            # Write out wavelength data
            fp.write(np.int64(len(data.wl)))
            fp.write(data.wl)

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_end = ptr + np.argmin(
                        np.logical_and(data.trial[ptr:] == trials[i],
                                       data.time[ptr:] == times[j]))

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.spec[i,:])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Convert wavelength data to FITS columns and make an HDU
            # from it; complication: astropy expects the dimensions of
            # the array to be (n_entries, n_wavelengths)
            nl = data.wl.shape[0]
            fmtstring = str(nl)+"D"
            wlcols = [fits.Column(name="Wavelength",
                                  format=fmtstring,
                                  unit="Angstrom", 
                                  array=data.wl.reshape(1,nl))]
            wlfits = fits.ColDefs(wlcols)
            wlhdu = fits.BinTableHDU.from_columns(wlcols)

            # Convert spectra to FITS columns, and make an HDU from
            # them
            speccols = []
            speccols.append(fits.Column(name="Trial", format="1K",
                                        unit="", array=data.trial))
            speccols.append(fits.Column(name="UniqueID", format="1K",
                                        unit="", array=data.id))
            speccols.append(fits.Column(name="Time", format="1D",
                                        unit="yr", array=data.time))
            speccols.append(fits.Column(name="L_lambda",
                                        format=fmtstring,
                                        unit="erg/s/A",
                                        array=data.spec))
            specfits = fits.ColDefs(speccols)
            spechdu = fits.BinTableHDU.from_columns(specfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
            hdulist.writeto(model_name+'_cluster_spec.fits', 
                            clobber=True)

    ################################################################
    # Write photometry file if we have the data for it
    ################################################################
    if 'phot' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_phot.txt', 'w')

            # Write header lines
            fp.write(("{:<18s}"*2).format('UniqueID', 'Time'))
            for f in data.filter_names:
                fp.write("{:<18s}".format(f))
            fp.write("\n")
            fp.write(("{:<18s}"*2).format('', '(yr)'))
            for f in data.filter_units:
                fp.write("({:s}".format(f)+")"+" "*(16-len(f)))
            fp.write("\n")
            nf = len(data.filter_names)
            fp.write(("{:<18s}"*2).
                     format('---------------', '---------------'))
            for i in range(nf):
                fp.write("{:<18s}".format('---------------'))
            fp.write("\n")

            # Write data
            for i in range(len(data.id)):
                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*((2+nf)*18-3)+"\n")
                fp.write("    {:11d}       {:11.5e}"
                         .format(data.id[i], data.time[i]))
                for j in range(nf):
                    fp.write("       {:11.5e}".format(data.phot[i,j]))
                fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_phot.bin', 'wb')

            # Write number of filters and filter names as ASCII
            nf = len(data.filter_names)
            fp.write(str(nf)+"\n")
            for i in range(nf):
                fp.write(data.filter_names[i] + " " + 
                         data.filter_units[i] + "\n")

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_end = ptr + np.argmin(
                        np.logical_and(data.trial[ptr:] == trials[i],
                                       data.time[ptr:] == times[j]))

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.phot[k,:])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
            cols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
            for i in range(len(data.filter_names)):
                cols.append(fits.Column(name=data.filter_names[i],
                                        unit=data.filter_units[i],
                                        format="1D",
                                        array=data.phot[:,i]))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_phot.fits',
                            clobber=True)
