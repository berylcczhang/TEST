"""
Function to write a set of output integrated files in SLUG2 format,
starting from an integrated data set as returned by
read_integrated. This can be used to translate data from one format to
another (e.g., bin to fits), or to consolidate multiple runs into a
single output file.
"""

import numpy as np
from cloudy import write_integrated_cloudylines
from cloudy import write_integrated_cloudyspec
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")


def write_integrated(data, model_name, fmt):
    """
    Function to write a set of output integrated files in SLUG2 format,
    starting from an integrated data set as returned by
    read_integrated.

    Parameters
    ----------
    data : namedtuple
       Integrated data to be written, in the namedtuple format returned
       by read_integrated
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
    if 'target_mass' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_prop.txt', 'w')

            # Write header lines
            fp.write(("{:<14s}"*8).
                     format('Time', 'TargetMass', 'ActualMass',
                            'LiveMass', 'ClusterMas', 'NumClusters',
                            'NumDisClust', 'NumFldStar') + "\n")
            fp.write(("{:<14s}"*8).
                     format('(yr)', '(Msun)', '(Msun)', '(Msun)',
                            '(Msun)', '', '', '') + "\n")
            fp.write(("{:<14s}"*8).
                     format('-----------', '-----------', '-----------',
                            '-----------', '-----------',
                            '-----------', '-----------',
                            '-----------') + "\n")

            # Write data
            ntime = len(data.time)
            ntrial = data.actual_mass.shape[-1]
            for i in range(ntrial):
                if i != 0:
                    fp.write("-"*(8*14-3)+"\n")
                for j in range(ntime):
                    fp.write(("{:11.5e}   {:11.5e}   {:11.5e}   " +
                              "{:11.5e}   {:11.5e}   {:11d}   " +
                              "{:11d}   {:11d}\n")
                             .format(data.time[j], 
                                     data.target_mass[j],
                                     data.actual_mass[j,i],
                                     data.live_mass[j,i],
                                     data.cluster_mass[j,i],
                                     data.num_clusters[j,i],
                                     data.num_dis_clusters[j,i],
                                     data.num_fld_stars[j,i]))

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_prop.bin', 'wb')

            # Write data
            ntime = len(data.time)
            ntrial = data.actual_mass.shape[-1]
            for i in range(ntrial):
                for j in range(ntime):
                    fp.write(data.time[j])
                    fp.write(data.target_mass[j])
                    fp.write(data.actual_mass[j,i])
                    fp.write(data.live_mass[j,i])
                    fp.write(data.cluster_mass[j,i])
                    fp.write(data.num_clusters[j,i])
                    fp.write(data.num_dis_clusters[j,i])
                    fp.write(data.num_fld_stars[j,i])

            # Close
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Figure out number of trials, and tile arrays
            ntrial = data.actual_mass.shape[-1]
            ntimes = len(data.time)
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            times = np.tile(data.time, ntrial)
            target_mass = np.tile(data.target_mass, ntrial)

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
            cols.append(fits.Column(name="TargetMass", format="1D",
                                    unit="Msun", array=target_mass))
            cols.append(
                fits.Column(name="ActualMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.actual_mass).flatten()))
            cols.append(
                fits.Column(name="LiveMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.live_mass).flatten()))
            cols.append(
                fits.Column(name="ClusterMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.cluster_mass).flatten()))
            cols.append(
                fits.Column(name="NumClusters", format="1K",
                            unit="", 
                            array=np.transpose(data.num_clusters).flatten()))
            cols.append(
                fits.Column(name="NumDisClust", format="1K",
                            unit="", 
                            array=np.transpose(data.num_dis_clusters).
                            flatten()))
            cols.append(
                fits.Column(name="NumFldStar", format="1K",
                            unit="", 
                            array=np.transpose(data.num_fld_stars).flatten()))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_integrated_prop.fits',
                            clobber=True)

    ################################################################
    # Write the spectra file if we have the data for it
    ################################################################
    if 'spec' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_spec.txt', 'w')

            # Write header lines
            fp.write(("{:<14s}"*3).
                     format('Time', 'Wavelength', 'L_lambda') + "\n")
            fp.write(("{:<14s}"*3).
                     format('(yr)', '(Angstrom)', '(erg/s/A)') + "\n")
            fp.write(("{:<14s}"*3).
                     format('-----------', '-----------', '-----------')
                     + "\n")

            # Write data
            ntime = len(data.time)
            ntrial = data.spec.shape[-1]
            nl = len(data.wl)
            for i in range(ntrial):
                if i != 0:
                    fp.write("-"*(3*14-3)+"\n")
                for j in range(ntime):
                    for k in range(nl):
                        fp.write("{:11.5e}   {:11.5e}   {:11.5e}\n"
                                 .format(data.time[j], data.wl[k],
                                         data.spec[k,j,i]))

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_spec.bin', 'wb')

            # Write out wavelength data
            fp.write(np.int64(len(data.wl)))
            fp.write(data.wl)

            # Write out times and spectra
            ntime = len(data.time)
            ntrial = data.spec.shape[-1]
            for i in range(ntrial):
                for j in range(ntime):
                    fp.write(data.time[j])
                    # This next line is needed to put the data into a
                    # contiguous block before writing
                    tmp = np.copy(data.spec[:,j,i])
                    fp.write(tmp)

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

            # Figure out number of trials, and tile arrays
            ntrial = data.actual_mass.shape[-1]
            ntimes = len(data.time)
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            times = np.tile(data.time, ntrial)

            # Convert spectra to FITS columns, and make an HDU from
            # them
            speccols = []
            speccols.append(fits.Column(name="Trial", format="1K",
                                        unit="", array=trial))
            speccols.append(fits.Column(name="Time", format="1D",
                                        unit="yr", array=times))
            speccols.append(fits.Column(name="L_lambda",
                                        format=fmtstring,
                                        unit="erg/s/A",
                                        array=np.transpose(data.spec).
                                        reshape(ntimes*ntrial, nl)))
            specfits = fits.ColDefs(speccols)
            spechdu = fits.BinTableHDU.from_columns(specfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
            hdulist.writeto(model_name+'_integrated_spec.fits', 
                            clobber=True)
                
    ################################################################
    # Write photometry file if we have the data for it
    ################################################################
    if 'phot' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_phot.txt', 'w')

            # Write header lines
            fp.write("{:<18s}".format('Time'))
            for f in data.filter_names:
                fp.write("{:<18s}".format(f))
            fp.write("\n")
            fp.write("{:<18s}".format('(yr)'))
            for f in data.filter_units:
                fp.write("({:s}".format(f)+")"+" "*(16-len(f)))
            fp.write("\n")
            fp.write("{:<18s}".format('---------------'))
            nf = len(data.filter_names)
            for i in range(nf):
                fp.write("{:<18s}".format('---------------'))
            fp.write("\n")

            # Write data
            ntime = data.phot.shape[1]
            ntrial = data.phot.shape[2]
            for i in range(ntrial):
                # Write separator between trials
                if i != 0:
                    fp.write("-"*((1+nf)*18-3)+"\n")
                for j in range(ntime):
                    fp.write("    {:11.5e}".format(data.time[j]))
                    for k in range(nf):
                        fp.write("       {:11.5e}".format(data.phot[k,j,i]))
                    fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_phot.bin', 'wb')

            # Write number of filters and filter names as ASCII
            nf = len(data.filter_names)
            fp.write(str(nf)+"\n")
            for i in range(nf):
                fp.write(data.filter_names[i] + " " + 
                         data.filter_units[i] + "\n")

            # Write data
            ntime = data.phot.shape[1]
            ntrial = data.phot.shape[2]
            for i in range(ntrial):
                for j in range(ntime):
                    fp.write(data.time[j])
                    # This next line is needed to put the data into a
                    # contiguous block before writing
                    tmp = np.copy(data.phot[:,j,i])
                    fp.write(tmp)

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Figure out number of trials, and tile arrays
            ntimes = data.phot.shape[1]
            ntrial = data.phot.shape[2]
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            times = np.tile(data.time, ntrial)
            nf = len(data.filter_names)

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
            for i in range(len(data.filter_names)):
                cols.append(
                    fits.Column(name=data.filter_names[i],
                                unit=data.filter_units[i],
                                format="1D",
                                array=np.transpose(data.phot[i,:,:]).
                                flatten()))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_integrated_phot.fits',
                            clobber=True)

    ################################################################
    # Write cloudy files if we have the data for them
    ################################################################
    if 'cloudy_inc' in data._fields:
        write_integrated_cloudyspec(data, model_name, fmt=fmt)
    if 'cloudy_linelum' in data._fields:
        write_integrated_cloudylines(data, model_name, fmt=fmt)
