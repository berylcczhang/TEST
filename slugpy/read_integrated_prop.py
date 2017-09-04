"""
Function to read a SLUG2 integrated_prop file.
"""

import numpy as np
from collections import namedtuple
from collections import defaultdict
import struct
import re
from .slug_open import slug_open

def read_integrated_prop(model_name, output_dir=None, fmt=None, 
                         verbose=False, read_info=None,
                         no_stellar_mass=False):
    """
    Function to read a SLUG2 integrated_prop file.

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
       no_stellar_mass : bool
          Prior to 7/15, output files did not contain the stellar_mass
          field; this can be detected automatically for ASCII and FITS
          formats, but not for binary format; if True, this specifies
          that the binary file being read does not contain a
          stellar_mass field; it has no effect for ASCII or FITS files

    Returns
       A namedtuple containing the following fields:

       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       target_mass : array, shape (N_times, N_trials)
          Target stellar mass at each time
       actual_mass : array, shape (N_times, N_trials)
          Actual mass of stars created up to each time in each trial
       live_mass : array, shape (N_times, N_trials)
          Mass of currently-alive stars at each time in each trial
       stellar_mass : array
          mass of all stars, living and stellar remnants
       cluster_mass : array, shape (N_times, N_trials)
          Stellar mass in non-disrupted clusters at each time in each
          trial
       num_clusters : array, shape (N_times, N_trials), dtype ulonglong
          Number of non-disrupted clusters present at each time in each
          trial
       num_dis_clusters : array, shape (N_times, N_trials), dtype ulonglong
          Number of disrupted clusters present at each time in each trial
       num_fld_stars : array, shape (N_times, N_trials), dtype ulonglong
          Number of living field stars (excluding those in disrupted 
          clusters and those being treated non-stochastically) present at
          each time in each trial
    """

    # Open file
    fp, fname = slug_open(model_name+"_integrated_prop", 
                          output_dir=output_dir,
                          fmt=fmt)

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False
        
    # Print status and record
    if verbose:
        print("Reading integrated properties for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare lists to hold data
    time = []
    target_mass = []
    actual_mass = []
    live_mass = []
    stellar_mass = []
    cluster_mass = []
    num_clusters = []
    num_dis_clusters = []
    num_fld_stars = []


    imf_is_var = False
    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()

        # Read the first header line
        hdr = fp.readline()

        # See if we have the stellar mass field; this was added later,
        # so we check in order to maintain backwards compatibility
        hdrsplit = hdr.split()
        if 'StellarMass' in hdrsplit:
            has_st_mass = True
        else:
            has_st_mass = False

        # Check for variable imf parameters
        if 'VP0' not in hdrsplit:
            imf_is_var = False
            checking_for_var = False
        else:
            imf_is_var = True
            checking_for_var = True
            vplist=[]
            p = 0;
            while (checking_for_var == True):
                
                if 'VP'+repr(p) in hdrsplit:
                    vplist.append(p)
                    p=p+1
                    vp_dict = defaultdict(list)

                else:
                  checking_for_var = False





        # Burn two header lines
        fp.readline()
        fp.readline()

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
            target_mass.append(float(data[1]))
            actual_mass.append(float(data[2]))
            live_mass.append(float(data[3]))
            if has_st_mass:
                stellar_mass.append(float(data[4]))
                cluster_mass.append(float(data[5]))
                num_clusters.append(int(data[6]))
                num_dis_clusters.append(int(data[7]))
                num_fld_stars.append(int(data[8]))
                datanumber=8
                
            else:
                stellar_mass.append(np.nan)
                cluster_mass.append(float(data[4]))
                num_clusters.append(int(data[5]))
                num_dis_clusters.append(int(data[6]))
                num_fld_stars.append(int(data[7]))
                datanumber=7

            if imf_is_var:
                for i in vplist:
                    datanumber += 1
                    vp_dict['VP'+repr(i)].append(float(data[datanumber]))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))
        
        # Read in number of variable parameters    
        firstdata = fp.read(struct.calcsize('i'))
        nvps = struct.unpack('i',firstdata)[0]
        
        if nvps > 0:
            imf_is_var = True
            vp_dict = defaultdict(list)

        # Suck file into memory
        data = fp.read()

        # Interpret data
        if not no_stellar_mass:
            datstr = 'LddddddQQQ'
            nfield = 10
        else:
            datstr = 'LdddddQQQ'
            nfield = 9

        if imf_is_var:
            datstr = datstr + nvps*'d'        
            nfield = len(datstr)
        
            
        nentry = len(data)//struct.calcsize(datstr)
        data_list = struct.unpack(datstr*nentry, data)

        # Stick data into correctly-named lists
        trial = data_list[0::nfield]
        time = data_list[1::nfield]
        target_mass = data_list[2::nfield]
        actual_mass = data_list[3::nfield]
        live_mass = data_list[4::nfield]
        if not no_stellar_mass:
            stellar_mass = data_list[5::nfield]
            cluster_mass = data_list[6::nfield]
            num_clusters = data_list[7::nfield]
            num_dis_clusters = data_list[8::nfield]
            num_fld_stars = data_list[9::nfield]
            lastindex=9
        else:
            cluster_mass = data_list[5::nfield]
            num_clusters = data_list[6::nfield]
            num_dis_clusters = data_list[7::nfield]
            num_fld_stars = data_list[8::nfield]
            stellar_mass = np.copy(np.array(cluster_mass))
            stellar_mass[...] = np.nan
            lastindex=8

        if imf_is_var:                  
            currentvp = 0
            while currentvp < nvps:  
                          
                vp_dict['VP'+repr(currentvp)].extend(data_list[lastindex+currentvp+1::nfield])
                currentvp += 1            
        


    elif fname.endswith('fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        target_mass = fp[1].data.field('TargetMass')
        actual_mass = fp[1].data.field('ActualMass')
        live_mass = fp[1].data.field('LiveMass')
        try:
            stellar_mass = fp[1].data.field('StellarMass')
        except KeyError:
            stellar_mass = np.zeros(live_mass.shape)
            stellar_mass[...] = np.nan
        cluster_mass = fp[1].data.field('ClusterMass')
        num_clusters = fp[1].data.field('NumClusters')
        num_dis_clusters = fp[1].data.field('NumDisClust')
        num_fld_stars = fp[1].data.field('NumFldStar')



        # Check for variable imf parameters
        imf_is_var = True
        if 'VP0' not in fp[1].data.columns.names:
            imf_is_var = False
            checking_for_var = False
        else:
            checking_for_var = True
            vplist=[]
            p = 0
            vp_dict = defaultdict(list)
            while (checking_for_var == True):
                if 'VP'+repr(p) in fp[1].data.columns.names:
                    
                    vp_dict['VP'+repr(p)].append(fp[1].data.field('VP'+repr(p)))         

                    p=p+1

                else:
                    checking_for_var = False



    # Close file
    fp.close()

    # Convert lists to arrays
    trial = np.array(trial)
    time = np.array(time)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    stellar_mass = np.array(stellar_mass)
    cluster_mass = np.array(cluster_mass)
    num_clusters = np.array(num_clusters, dtype='ulonglong')
    num_dis_clusters = np.array(num_dis_clusters, dtype='ulonglong')
    num_fld_stars = np.array(num_fld_stars, dtype='ulonglong')
    if imf_is_var:               
        for VPn in vp_dict:         
            vp_dict[VPn]=np.array(vp_dict[VPn])
            if fname.endswith('.fits'):
                vp_dict[VPn]=vp_dict[VPn][0]
                


    # Figure out if we have a number of trials with identical times,
    # indicating fixed output times, or if each trial has random times;
    # reshape time array appropriately
    ntrial = len(np.unique(trial))
    ntime = len(time)//ntrial
    if ntime != len(time):
        if np.amin(time[:ntime] == time[ntime:2*ntime]):
            time = time[:ntime]

    # Prune / reshape the output arrays
    target_mass = np.transpose(target_mass.reshape(ntrial, ntime))
    actual_mass = np.transpose(actual_mass.reshape(ntrial, ntime))
    live_mass = np.transpose(live_mass.reshape(ntrial, ntime))
    stellar_mass = np.transpose(stellar_mass.reshape(ntrial, ntime))
    cluster_mass = np.transpose(cluster_mass.reshape(ntrial, ntime))
    num_clusters = np.transpose(num_clusters.reshape(ntrial, ntime))
    num_dis_clusters = np.transpose(num_dis_clusters.reshape(ntrial, ntime))
    num_fld_stars = np.transpose(num_fld_stars.reshape(ntrial, ntime))


    if imf_is_var:               
        for VPn in vp_dict:         
            vp_dict[VPn]=np.transpose(vp_dict[VPn].reshape(ntrial,ntime))




    # Build the namedtuple to hold output
    out_type = namedtuple('integrated_prop',
                          ['time', 'target_mass', 'actual_mass',
                           'live_mass', 'stellar_mass',
                           'cluster_mass', 'num_clusters',
                           'num_dis_clusters', 'num_fld_stars'])
    out = out_type(time, target_mass, actual_mass, live_mass,
                   stellar_mass, cluster_mass, num_clusters, 
                   num_dis_clusters, num_fld_stars)
                   
    if imf_is_var:
        vpn_tuple = ()
        
        for VPn in vp_dict:         

            out_type = namedtuple('cluster_prop',out_type._fields+(VPn,))
            vpn_tuple = vpn_tuple + (vp_dict[VPn],)
                                
            out = out_type(time, target_mass, actual_mass, live_mass,
                   stellar_mass, cluster_mass, num_clusters, 
                   num_dis_clusters, num_fld_stars,*vpn_tuple)                     

    # Return
    return out

