"""
Function to read a SLUG2 cluster_yield file.
"""

import numpy as np
from collections import namedtuple
import struct
from slug_open import slug_open

def read_cluster_yield(model_name, output_dir=None, fmt=None, 
                      verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_spec file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the file to be read. If one of these is set, the
          function will only attempt to open ASCII-('txt' or 'ascii'), 
          binary ('bin' or 'binary'), or FITS ('fits' or 'fits2')
          formatted output, ending in .txt., .bin, or .fits,
          respectively. If set to None, the code will try to open
          ASCII files first, then if it fails try binary files, and if
          it fails again try FITS files.
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       name : array of strings, shape (N_iso)
          Atomic symbols of the isotopes included in the yield table
       Z : array of int, shape (N_iso)
          Atomic numbers of the isotopes included in the yield table
       A : array of int, shape (N_iso)
          Atomic mass number of the isotopes included in the yield table
       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          times at which cluster spectra are output, in yr
       yld : array, shape (N_cluster, N_iso)
          Yield of each isotope, defined as the instantaneous amount
          produced up to that time; for unstable isotopes, this
          includes the effects of decay since production
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_yield", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster yield for model "+model_name)
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

        # Prepare storage
        isotope_name = []
        isotope_Z = []
        isotope_A = []
        yldtmp = []
        trialptr = 0
        id_save = -1
        time_save = -1.0

        # Read the first cluster
        for entry in fp:

            # Have we hit a separator?
            if entry[:3] == '---':
                trialptr += 1
                break

            # Parse this line
            data = entry.split()
            cl_id = long(data[0])
            t = float(data[1])
            iso_name = data[2].title()
            iso_Z = int(data[3])
            iso_A = int(data[4])
            iso_yld = float(data[5])

            # Are we still on the same cluster? If not, break
            if id_save < 0:
                id_save = cl_id
                time_save = t
            elif (id_save != cl_id) or (t != time_save):
                break

            # If we're here, we're still on the same cluster, so save
            # the isotope information
            isotope_name.append(iso_name)
            isotope_Z.append(iso_Z)
            isotope_A.append(iso_A)
            yldtmp.append(iso_yld)

        # If we're here, we're read the first block of
        # isotopes. Append the time, cluster ID, trial, and yield to
        # the list
        cluster_id = [id_save]
        time = [time_save]
        trial = [0]
        yld = [yldtmp]

        # Prepare for loop over rest of file
        if trialptr == 0:
            yldtmp = [iso_yld]
            isoptr = 1
        else:
            yldtmp = []
            isoptr = 0

        # Loop through rest of file
        for entry in fp:

            # Save if isotope pointer has looped
            if isoptr == len(isotope_name):
                cluster_id.append(cl_id)
                time.append(t)
                yld.append(yldtmp)
                trial.append(trialptr)
                yldtmp = []
                isoptr = 0

            # Have we hit a separator?
            if entry[:3] == '---':
                trialptr += 1
                continue

            # Parse this line
            data = entry.split()
            cl_id = long(data[0])
            t = float(data[1])
            yldtmp.append(float(data[5]))

            # Increment isotope pointer
            isoptr = isoptr + 1

        # Save last block
        cluster_id.append(cl_id)
        time.append(t)
        yld.append(yldtmp)
        trial.append(trialptr)
        
        # Convert to arrays
        isotope_name = np.array(isotope_name)
        isotope_Z = np.array(isotope_Z, dtype=int)
        isotope_A = np.array(isotope_A, dtype=int)
        time = np.array(time)
        cluster_id = np.array(cluster_id)
        yld = np.array(yld)
        trial = np.array(trial, dtype=int)

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read number of isotopes
        niso = struct.unpack('L', fp.read(struct.calcsize('L')))[0]

        # Read isotope data
        data = fp.read(struct.calcsize(('c'*4+'II')*niso))
        data_list = struct.unpack(('c'*4+'II')*niso, data)
        isotope_name = np.array(
            [ (data_list[6*i]+data_list[6*i+1]+
               data_list[6*i+2]+data_list[6*i+3]).strip().
              title() for i in range(niso) ])
        isotope_Z = np.array(data_list[4::6], dtype=int)
        isotope_A = np.array(data_list[5::6], dtype=int)

        # Now read remainder of file
        buf = fp.read()

        # Prepare storage
        time = []
        trial = []
        cluster_id = np.zeros(0, dtype=int)
        yld = np.zeros((0,niso))

        # Parse the buffer
        ptr = 0
        hdrstr = 'LdL'
        hdrsize = struct.calcsize('LdL')
        while ptr < len(buf):

            # Read number of clusters in this block
            trialptr, t, ncluster \
                = struct.unpack(hdrstr, buf[ptr:ptr+hdrsize])

            # Skip if no clusters
            if ncluster == 0:
                continue

            # Add to time and trial arrays
            time = time + [t]*ncluster
            trial = trial + [trialptr]*ncluster

            # Read ID's and yields for these clusters
            blockstr = ('L'+'d'*niso)*ncluster
            blocksize = struct.calcsize(blockstr)
            block_data = struct.unpack(
                blockstr, buf[ptr+hdrsize:ptr+hdrsize+blocksize])
            cl_id = np.array(block_data[::niso+1], dtype=int)
            yldtmp = np.array([block_data[(niso+1)*i+1:(niso+1)*(i+1)]
                               for i in range(ncluster)])

            # Add to arrays
            cluster_id = np.append(cluster_id, cl_id)
            yld = np.append(yld, yldtmp, axis=0)

            # Move pointer
            ptr = ptr + hdrsize + blocksize

        # Convert to arrays
        time = np.array(time)
        trial = np.array(trial)

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Read data
        isotope_name = fp[1].data['Name']
        isotope_Z = fp[1].data['Z']
        isotope_A = fp[1].data['A']
        cluster_id = fp[2].data['UniqueID']
        trial = fp[2].data['Trial']
        time = fp[2].data['Time']
        yld = fp[2].data['Yield']

    # Close file
    fp.close()

    # Build output holder
    fieldnames = ['name', 'Z', 'A', 'cluster_id', 'time', 'trial', 'yld']
    fields = [ isotope_name, isotope_Z, isotope_A, cluster_id, time,
               trial, yld]
    out_type = namedtuple('integrated_yield', fieldnames)
    out = out_type(*fields)

    # Return
    return out

