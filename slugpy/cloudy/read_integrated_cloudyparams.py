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
       cloudy_r1 : array
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
          simulated by cloudy
       cloudy_U0 : array
          ionization parameter of the HII reegion at its inner edge
       cloudy_Omega : array
          Yeh & Matzner (2012) wind parameter for the HII region
          simulated by cloudy
       cloudy_zeta : array
          Krumholz & Matzner (2009) radiation pressure parameter for
          the HII region, again approximate; values of zeta >~1
          indicate that radiation pressure is dynamically important

       The following fields may or may not be present, depending on
       what is found in the output file:
       
       cloudy_hden_out : array
          volume-averaged number density produced by the cloudy
          calculation
       cloudy_r1_out : array
          HII region outer radius produced by cloudy
       cloudy_Omega_out : array
          value of Omega computed using cloudy_r1_out
       cloudy_zeta_out : array
          value of zeta computed using cloudy_r1_out

    Notes

       The relationships between U, U0, Omega, r0, r1, hden, and
       QH0 used in deriving various parameters are valid only in
       the limit of negligible radiation pressure. They may be
       significantly off if radiation pressure is significant,
       i.e., if zeta >~ 1.
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

        # Read the first header line and see what fields are present
        line = fp.readline()
        if 'Hden_out' in line: hden_out_set = True
        else: hden_out_set = False
        if 'R1_out' in line: r1_out_set = True
        else: r1_out_set = False
        if 'Omega_out' in line: Omega_out_set = True
        else: Omega_out_set = False
        if 'Zeta_out' in line: zeta_out_set = True
        else: zeta_out_set = False
                
        # Burn two header lines
        line = fp.readline()
        line = fp.readline()

        # Prepare output holders
        time = []
        trial = []
        hden = []
        r0 = []
        r1 = []
        QH0 = []
        covFac = []
        U = []
        U0 = []
        Omega = []
        zeta = []
        if hden_out_set: hden_out = []
        if r1_out_set: r1_out = []
        if Omega_out_set: Omega_out = []
        if zeta_out_set: zeta_out = []

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
            r1.append(float(data[3]))
            QH0.append(float(data[4]))
            covFac.append(float(data[5]))
            U.append(float(data[6]))
            U0.append(float(data[7]))
            Omega.append(float(data[8]))
            zeta.append(float(data[9]))
            ptr = 10
            if hden_out_set:
                hden_out.append(float(data[ptr]))
                ptr += 1
            if r1_out_set:
                r1_out.append(float(data[ptr]))
                ptr += 1
            if Omega_out_set:
                Omega_out.append(float(data[ptr]))
                ptr += 1
            if zeta_out_set:
                zeta_out.append(float(data[ptr]))
                ptr += 1

        # Convert to arrays
        time = np.array(time)
        trial = np.array(trial)
        hden = np.array(hden)
        r0 = np.array(r0)
        r1 = np.array(r1)
        QH0 = np.array(QH0)
        covFac = np.array(covFac)
        U = np.array(U)
        U0 = np.array(U0)
        Omega = np.array(Omega)
        zeta = np.array(zeta)
        if hden_out_set: hden_out = np.array(hden_out, dtype='float')
        if r1_out_set: r1_out = np.array(r1_out, dtype='float')
        if Omega_out_set: Omega_out = np.array(Omega_out, dtype='float')
        if zeta_out_set: zeta_out = np.array(zeta_out, dtype='float')

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read leading bits to see which fields we have
        data = fp.read(struct.calcsize('bbbb'))
        field_bits = struct.unpack('bbbb', data)
        hden_out_set = field_bits[0] != 0
        r1_out_set = field_bits[1] != 0
        Omega_out_set = field_bits[2] != 0
        zeta_out_set = field_bits[3] != 0

        # Suck rest of file into memory
        data = fp.read()

        # Break data up into entries
        datastr = 'Ldddddddddd'
        if hden_out_set: datastr += 'd'
        if r1_out_set: datastr += 'd'
        if Omega_out_set: datastr += 'd'
        if zeta_out_set: datastr += 'd'
        nfield = len(datastr)
        nentry = len(data)/struct.calcsize(datastr)
        data_list = struct.unpack(datastr*nentry, data)

        # Put data into arrays
        trial = np.array(data_list[0::nfield])
        time = np.array(data_list[1::nfield])
        hden = np.array(data_list[2::nfield])
        r0 = np.array(data_list[3::nfield])
        r1 = np.array(data_list[4::nfield])
        QH0 = np.array(data_list[5::nfield])
        covFac = np.array(data_list[6::nfield])
        U = np.array(data_list[7::nfield])
        U0 = np.array(data_list[8::nfield])
        Omega = np.array(data_list[9::nfield])
        zeta = np.array(data_list[10::nfield])
        ptr = 11
        if hden_out_set:
            hden_out.extend(data_list[ptr::nfield])
            ptr += 1
        if r1_out_set:
            r1_out.extend(data_list[ptr::nfield])
            ptr += 1
        if Omega_out_set:
            Omega_out.extend(data_list[ptr::nfield])
            ptr += 1
        if zeta_out_set:
            zeta_out.extend(data_list[ptr::nfield])
            ptr += 1
        
    elif fmt == 'fits' or fname.endswith('fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        hden = fp[1].data.field('Hden')
        r0 = fp[1].data.field('R0')
        r1 = fp[1].data.field('R1')
        QH0 = fp[1].data.field('QH0')
        covFac = fp[1].data.field('covFac')
        U = fp[1].data.field('U')
        U0 = fp[1].data.field('U0')
        Omega = fp[1].data.field('Omega')
        zeta = fp[1].data.field('zeta')
        try:
            hden_out = fp[1].data.field('Hden_out')
            hden_out_set = True
        except KeyError:
            hden_out_set = False
        try:
            r1_out = fp[1].data.field('R1_out')
            r1_out_set = True
        except KeyError:
            r1_out_set = False
        try:
            Omega_out = fp[1].data.field('Omega_out')
            Omega_out_set = True
        except KeyError:
            Omega_out_set = False
        try:
            zeta_out = fp[1].data.field('zeta_out')
            zeta_out_set = True
        except KeyError:
            zeta_out_set = False
        
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
    r1 = np.transpose(r1.reshape(ntrial, ntime))
    QH0 = np.transpose(QH0.reshape(ntrial, ntime))
    covFac = np.transpose(covFac.reshape(ntrial, ntime))
    U = np.transpose(U.reshape(ntrial, ntime))
    U0 = np.transpose(U0.reshape(ntrial, ntime))
    Omega = np.transpose(Omega.reshape(ntrial, ntime))
    zeta = np.transpose(zeta.reshape(ntrial, ntime))
    if hden_out_set:
        hden_out = np.transpose(hden_out.reshape(ntrial, ntime))
    if r1_out_set:
        hden_out = np.transpose(r1_out.reshape(ntrial, ntime))
    if Omega_out_set:
        Omega_out = np.transpose(Omega_out.reshape(ntrial, ntime))
    if zeta_out_set:
        zeta_out = np.transpose(zeta_out.reshape(ntrial, ntime))

    # Put output into namedtuple
    fields = ['time', 'cloudy_hden', 'cloudy_r0', 'cloudy_r1',
              'cloudy_QH0', 'cloudy_covFac', 'cloudy_U',
              'cloudy_U0', 'cloudy_Omega', 'cloudy_zeta']
    data = [time, hden, r0, r1, QH0, covFac, U, U0, Omega, zeta]
    if hden_out_set:
        fields.append('cloudy_hden_out')
        data.append(hden_out)
    if r1_out_set:
        fields.append('cloudy_r1_out')
        data.append(r1_out)
    if Omega_out_set:
        fields.append('cloudy_Omega_out')
        data.append(Omega_out)
    if zeta_out_set:
        fields.append('cloudy_zeta_out')
        data.append(zeta_out)
    cloudyparams_type = namedtuple('integrated_cloudyparams', fields)
    out = cloudyparams_type(*data)

    # Return
    return out
