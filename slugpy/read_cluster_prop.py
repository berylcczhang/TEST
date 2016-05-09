"""
Function to read a SLUG2 cluster_prop file.
"""

import numpy as np
from collections import namedtuple
from collections import defaultdict
import struct
from slug_open import slug_open

def read_cluster_prop(model_name, output_dir=None, fmt=None, 
                      verbose=False, read_info=None,
                      no_stellar_mass=False):
    """
    Function to read a SLUG2 cluster_prop file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the output is located; if set to None,
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

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          time at which cluster's properties are being evaluated
       form_time : array
          time when cluster formed
       lifetime : array
          time at which cluster will disrupt
       target_mass : array
          target cluster mass
       actual_mass : array
          actual mass at formation
       live_mass : array
          mass of currently living stars
       stellar_mass : array
          mass of all stars, living and stellar remnants
       num_star : array, dtype ulonglong
          number of living stars in cluster being treated stochastically
       max_star_mass : array
          mass of most massive living star in cluster
       A_V : array
          A_V value for each cluster, in mag (present only if SLUG was
          run with extinction enabled)
       vpn_tuple : contains arrays for any variable parameters we have (eg: VP0,
            VP1,VP2...) in the IMF.
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_prop",
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster properties for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Prepare lists to hold data
    cluster_id = []
    trial = []
    time = []
    form_time = []
    lifetime = []
    target_mass = []
    actual_mass = []
    live_mass = []
    stellar_mass = []
    num_star = []
    max_star_mass = []
    
    imf_is_var = False
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the first header line
        hdr = fp.readline()

        # See if we have extinction
        hdrsplit = hdr.split()
        if 'A_V' in hdrsplit:
            extinct = True
            A_V = []
        else:
            extinct = False

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
                
                if 'VP'+`p` in hdrsplit:
                    vplist.append(p)
                    p=p+1
                    vp_dict = defaultdict(list)

                else:
                  checking_for_var = False
                    
                
                   

        # See if we have the stellar mass field; this was added later,
        # so we check in order to maintain backwards compatibility
        if 'StellarMass' in hdrsplit:
            has_st_mass = True
        else:
            has_st_mass = False

        # Burn the next two header lines
        fp.readline()
        fp.readline()

        # Read data
        trialptr = 0
        for entry in fp:
            if entry[:3] == '---':
                # Separator line
                trialptr = trialptr + 1
                continue
            data = entry.split()
            cluster_id.append(long(data[0]))
            trial.append(trialptr)
            time.append(float(data[1]))
            form_time.append(float(data[2]))
            lifetime.append(float(data[3]))
            target_mass.append(float(data[4]))
            actual_mass.append(float(data[5]))
            live_mass.append(float(data[6]))
            if has_st_mass:
                stellar_mass.append(float(data[7]))
                num_star.append(long(data[8]))
                max_star_mass.append(float(data[9]))
                if extinct:
                    A_V.append(float(data[10]))
                   
                #Read in variable parameter values
                if extinct:
                    datanumber = 10
                else:
                    datanumber = 9
                if imf_is_var:
                    for i in vplist:
                        
                        datanumber += 1
                        vp_dict['VP'+`i`].append(float(data[datanumber]))

            else:
                stellar_mass.append(np.nan)
                num_star.append(long(data[7]))
                max_star_mass.append(float(data[8]))
                if extinct:
                    A_V.append(float(data[9]))
                #Read in variable parameter values
                if imf_is_var:
                    if extinct:
                        datanumber = 9
                    else:
                        datanumber = 8
                    for i in vplist:
                       
                        datanumber += 1
                        vp_dict['VP'+`i`].append(float(data[datanumber]))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read the first byte to see if we have extinction turned on
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0
        if extinct:
            A_V = []
            
        # Read the next few bytes to see if we have any variable parameters
        data = fp.read(struct.calcsize('i'))
        nvps = struct.unpack('i', data)[0]

        if nvps > 0:
            imf_is_var = True
            vp_dict = defaultdict(list)

        # Go through file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('LdL'))
            if len(data) < struct.calcsize('LdL'):
                break
            trialptr, t, ncluster = struct.unpack('LdL', data)

            # Skip if no clusters
            if ncluster == 0:
                continue

            # Read the next block of clusters
            if no_stellar_mass:
                datastr = 'LdddddQd'

            else:
                datastr = 'LddddddQd'
            if extinct:
                datastr = datastr+'d'
            
            if imf_is_var:
                datastr = datastr + nvps*'d'

            data = fp.read(struct.calcsize(datastr)*ncluster)
            data_list = struct.unpack(datastr*ncluster, data)

            # Pack these clusters into the data list
            if extinct:
                nfield = 10
            else:
                nfield = 9
            if no_stellar_mass:
                nfield = nfield-1
            if imf_is_var:
                nfield = len(datastr)            
            cluster_id.extend(data_list[0::nfield])
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)
            form_time.extend(data_list[1::nfield])
            lifetime.extend(data_list[2::nfield])
            target_mass.extend(data_list[3::nfield])
            actual_mass.extend(data_list[4::nfield])
            live_mass.extend(data_list[5::nfield])
            if no_stellar_mass:
                stellar_mass.extend([np.nan]*(len(data_list)/nfield))
                num_star.extend(data_list[5::nfield])
                max_star_mass.extend(data_list[7::nfield])
                if extinct:
                    A_V.extend(data_list[8::nfield])
                    if imf_is_var:                  
                        currentvp = 0
                        while currentvp < nvps:                
                            vp_dict['VP'+`currentvp`].extend(data_list[8+currentvp+1::nfield])
                            currentvp += 1
                elif imf_is_var:                  
                    currentvp = 0
                    while currentvp < nvps:                
                        vp_dict['VP'+`currentvp`].extend(data_list[7+currentvp+1::nfield])
                        currentvp += 1
                        
                    
            else:
                stellar_mass.extend(data_list[6::nfield])
                num_star.extend(data_list[7::nfield])
                max_star_mass.extend(data_list[8::nfield])
                if extinct:
                    A_V.extend(data_list[9::nfield])
                    if imf_is_var:
                        currentvp = 0
                        while currentvp < nvps:                
                            vp_dict['VP'+`currentvp`].extend(data_list[9+currentvp+1::nfield])
                            currentvp += 1      
                elif imf_is_var:                  
                    currentvp = 0
                    while currentvp < nvps:                
                        vp_dict['VP'+`currentvp`].extend(data_list[8+currentvp+1::nfield])
                        currentvp += 1                            
                              
    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        form_time = fp[1].data.field('FormTime')
        lifetime = fp[1].data.field('Lifetime')
        target_mass = fp[1].data.field('TargetMass')
        actual_mass = fp[1].data.field('BirthMass')
        live_mass = fp[1].data.field('LiveMass')
        try:
            stellar_mass = fp[1].data.field('StellarMass')
        except KeyError:
            stellar_mass = np.zeros(live_mass.shape)
            stellar_mass[...] = np.nan
        num_star = fp[1].data.field('NumStar')
        max_star_mass = fp[1].data.field('MaxStarMass')
        if 'A_V' in fp[1].data.columns.names:
            extinct = True
            A_V = fp[1].data.field('A_V')
        else:
            extinct = False
            
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
                if 'VP'+`p` in fp[1].data.columns.names:
                    
                    vp_dict['VP'+`p`].append(fp[1].data.field('VP'+`p`))         

                    p=p+1

                else:
                    checking_for_var = False
    # Close file
    fp.close()

    # Convert lists to arrays
    cluster_id = np.array(cluster_id, dtype='uint')
    trial = np.array(trial, dtype='uint')
    time = np.array(time)
    form_time = np.array(form_time)
    lifetime = np.array(lifetime)
    target_mass = np.array(target_mass)
    actual_mass = np.array(actual_mass)
    live_mass = np.array(live_mass)
    stellar_mass = np.array(stellar_mass)
    num_star = np.array(num_star, dtype='ulonglong')
    max_star_mass = np.array(max_star_mass)
    if extinct:
        A_V = np.array(A_V)
    if imf_is_var:        
        
        for VPn in vp_dict:         
            vp_dict[VPn]=np.array(vp_dict[VPn])

            if fname.endswith('.fits'):
                vp_dict[VPn]=vp_dict[VPn][0]
                
            
        
        
    # Build the namedtuple to hold output
    if extinct:
        out_type = namedtuple('cluster_prop',
                              ['id', 'trial', 'time', 'form_time', 
                               'lifetime', 'target_mass', 'actual_mass', 
                               'live_mass', 'stellar_mass', 'num_star', 
                               'max_star_mass', 'A_V'])
        out = out_type(cluster_id, trial, time, form_time, lifetime, 
                       target_mass, actual_mass, live_mass, stellar_mass,
                       num_star, max_star_mass, A_V)
        
        if imf_is_var:
            vpn_tuple = ()

            for VPn in vp_dict:         
                out_type = namedtuple('cluster_prop',out_type._fields+(VPn,))
                vpn_tuple = vpn_tuple + (vp_dict[VPn],)

            out = out_type(cluster_id, trial, time, form_time, lifetime, 
                       target_mass, actual_mass, live_mass, stellar_mass,
                       num_star, max_star_mass, A_V,*vpn_tuple)    
            
                        
                  

                            
                       
    else:
        out_type = namedtuple('cluster_prop',
                              ['id', 'trial', 'time', 'form_time', 
                               'lifetime', 'target_mass', 'actual_mass', 
                               'live_mass', 'stellar_mass', 'num_star', 
                               'max_star_mass'])
        out = out_type(cluster_id, trial, time, form_time, lifetime, 
                       target_mass, actual_mass, live_mass, stellar_mass,
                       num_star, max_star_mass)
        if imf_is_var:
            vpn_tuple = ()

            for VPn in vp_dict:         

                out_type = namedtuple('cluster_prop',out_type._fields+(VPn,))
                vpn_tuple = vpn_tuple + (vp_dict[VPn],)

                                
            out = out_type(cluster_id, trial, time, form_time, lifetime, 
                       target_mass, actual_mass, live_mass, stellar_mass,
                       num_star, max_star_mass,*vpn_tuple)   
    # Return
    return out

