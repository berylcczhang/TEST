"""
Function to read a SLUG2 integrated_cloudyparams file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open

def read_integrated_cloudyparams(model_name, output_dir=None, fmt=None,
                                 verbose=False, read_info=None):
    """
    Function to read a SLUG2 integrated_cloudyparams file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : string
          Format for the file to be read. Allowed values are 'ascii',
          'bin' or 'binary, and 'fits'. If one of these is set, the code
          will only attempt to open ASCII-, binary-, or FITS-formatted
          output, ending in .txt., .bin, or .fits, respectively. If set
          to None, the code will try to open ASCII files first, then if
          it fails try binary files, and if it fails again try FITS
          files.
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       cloudy_hden : array
          number density of H nuclei at the inner edge of the ionized
          region simulated by cloudy
       cloudy_r0 : array
          inner radius of the ionized region simulated by cloudy
       cloudy_rS : array
          outer radius of the ionized region simulated by cloudy (approximate!)
       cloudy_QH0 : array
          ionizing luminosity used in the cloudy computation
       cloudy_covFac : array
          covering factor assumed in the cloudy computation; only a
          fraction covFac of the ionizing photons are assumed to
          produce emission within the HII region, while the remainder
          are assumed to escape
       cloudy_U : array
          volume-averaged ionization parameter of the HII region
          simulated by cloudy; note that this value is approximate,
          not exact, and the approximation can be very poor if
          radiation pressure effects are significant
       cloudy_Omega : array
          Yeh & Matzner (2012) wind parameter for the HII region
          simulated by cloudy; as with U, this value is approximate,
          and the approximation is valid only if radiation pressure
          effects are small
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_cloudyparams", 
                          output_dir=output_dir,
                          fmt=fmt)
    if read_info is not None:
        read_info['fname'] = fname

    # Print status
    if verbose:
        print("Reading integrated cloudy parameters for "
              "model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Burn 3 header lines
        hdr = fp.readline()
        hdr = fp.readline()
        hdr = fp.readline()

        # Prepare output holders
        time = []
        trial = []
        hden = []
        r0 = []
        rS = []
        QH0 = []
        covFac = []
        U = []
        Omega = []

        # Read data
        trial = []
        trialptr = 0
        for entry in fp:
            if entry[:3] == '---':
                trialptr = trialptr+1
                continue       # Skip separator lines
            trial.append(trialptr)
            data = entry.split()
            time.append(float(data[0]))
            hden.append(float(data[1]))
            r0.append(float(data[2]))
            rS.append(float(data[3]))
            QH0.append(float(data[4]))
            covFac.append(float(data[5]))
            U.append(float(data[6]))
            Omega.append(float(data[7]))

        # Convert to arrays
        time = np.array(time)
        trial = np.array(trial)
        hden = np.array(hden)
        r0 = np.array(r0)
        rS = np.array(rS)
        QH0 = np.array(QH0)
        covFac = np.array(covFac)
        U = np.array(U)
        Omega = np.array(Omega)

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Suck file into memory
        data = fp.read()

        # Break data up into entries
        datastr = 'Ldddddddd'
        nfield = len(datastr)
        nentry = len(data)/struct.calcsize(datastr)
        data_list = struct.unpack(datastr*nentry, data)

        # Put data into arrays
        trial = np.array(data_list[0::nfield])
        time = np.array(data_list[1::nfield])
        hden = np.array(data_list[2::nfield])
        r0 = np.array(data_list[3::nfield])
        rS = np.array(data_list[4::nfield])
        QH0 = np.array(data_list[5::nfield])
        covFac = np.array(data_list[6::nfield])
        U = np.array(data_list[7::nfield])
        Omega = np.array(data_list[8::nfield])
        
    elif fmt == 'fits':

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        hden = fp[1].data.field('Hden')
        r0 = fp[1].data.field('R0')
        rS = fp[1].data.field('RS')
        QH0 = fp[1].data.field('QH0')
        covFac = fp[1].data.field('covFac')
        U = fp[1].data.field('U')
        Omega = fp[1].data.field('Omega')
        
    # Close file
    fp.close()

    # Figure out number of times and trials; if trials have identical
    # times, remove the duplicate information
    ntrial = len(np.unique(trial))
    ntime = len(time)/ntrial
    if ntime != len(time):
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]

    # Reshape the output arrays
    hden = np.transpose(hden.reshape(ntrial, ntime))
    r0 = np.transpose(r0.reshape(ntrial, ntime))
    rS = np.transpose(rS.reshape(ntrial, ntime))
    QH0 = np.transpose(QH0.reshape(ntrial, ntime))
    covFac = np.transpose(covFac.reshape(ntrial, ntime))
    U = np.transpose(U.reshape(ntrial, ntime))
    Omega = np.transpose(Omega.reshape(ntrial, ntime))
    
    # Put output into namedtuple
    cloudyparams_type \
        = namedtuple('integrated_cloudyparams',
                     ['time', 'cloudy_hden', 'cloudy_r0', 'cloudy_rS',
                      'cloudy_QH0', 'cloudy_covFac', 'cloudy_U',
                      'cloudy_Omega'])
    out = cloudyparams_type(time, hden, r0, rS, QH0, covFac, U, Omega)

    # Return
    return out
