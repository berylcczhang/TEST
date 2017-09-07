"""
This script takes outputs of several slug runs performed with
identical parameters and combines them into a single monolithic output
file. This capability is usually used to join together outputs that
have been run on different machines or different nodes. Functionally
it is similar to consolidate.py, but rather than assuming the naming
convention follows slug's defaults for parallel runs, it requires the
user to specify all the input and output file names.
"""

import argparse
import re
import os
import os.path as osp
import copy
import glob
import sys
import warnings
try:
    from slugpy import read_cluster, combine_cluster, \
        read_integrated, combine_integrated, write_cluster, \
        write_integrated
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import read_cluster, combine_cluster, \
            read_integrated, combine_integrated, write_cluster, \
            write_integrated
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Step 1: grab command line arguments
parser = argparse. \
         ArgumentParser(
             description=
             "Script to join multiple slug outputs into a single output file")
parser.add_argument("name", nargs='+',
                    help="base name(s) or file name(s) of runs "
                    "to be combined")
parser.add_argument("-o", "--outname", default="combined",
                    help="base name of the combined output file")
parser.add_argument("--delete", default=False, action="store_true",
                    help="delete original files after consolidation")
parser.add_argument("-f", "--fmt", default=None,
                    help="format for combined output; valid "+
                    "options are ascii, bin, fits, and fits2; "+
                    "default is same format as first input file")
parser.add_argument("--halt", default=False,
                    action="store_true",
                    help="halt execution on read failure; if this"
                    "option is not set, a warning is issued instead, "
                    "and processing continues")
parser.add_argument("-v", "--verbose", action="store_true",
                    default=False, help="produce verbose output")
args = parser.parse_args()

# Step 2: construct list of runs to be combined; the user can provide
# either file names or base names, so this is a bit tricky. We look
# for any files matching the specified names, and we also look for any
# files matching those names with an added extension corresponding to
# a slug output file. Once we have the complete list of file names, we
# identify ones that follow the name format for a slug output file and
# construct a set of unique base names.
extensions = ["_cluster_*.txt", "_cluster_*.bin", "_cluster_*.fits",
              "_integrated_*.txt", "_integrated_*.bin",
              "_integrated_*.fits"]
candidates = []
for n in args.name:
    # Add files matching input names to candidate file list
    candidates = candidates + glob.glob(n)
    # Add extensions to candidate file list
    for e in extensions:
        candidates = candidates + glob.glob(n+e)
basenames = []
for c in candidates:
    # Search for valid file names
    for e in extensions:
        e_re = e.replace(".", "\.").replace("*", ".*") # regex syntax
        if re.search(e_re, c) is not None:
            bname = re.sub(e_re, "", c)
            if bname not in basenames:
                basenames.append(bname)
if args.verbose:
    print("Combining the following runs:")
    for b in basenames:
        print(b)

# Step 3: read all cluster data into memory
fmt = args.fmt
files_read = []
data = []
for b in basenames:
    cl_files = glob.glob(b + "_cluster_*")
    if len(cl_files) > 0:
        rinfo={}
        try:
            if args.verbose:
                print("Reading cluster data from run "+b)
            cldata = read_cluster(b, nofilterdata=True, read_info=rinfo)
            data.append(cldata)
            for k in rinfo.keys():
                if k != 'format':
                    files_read.append(rinfo[k])
                elif fmt is None:
                    fmt = rinfo['format']
        except IOError as e:
            if args.halt:
                raise IOError(e)
            else:
                warnings.warn("IOError during processing of cluster "
                              "data for run "+ b + ":" 
                              + str(e) + "; processing continues", 
                              RuntimeWarning)

# Step 4: combine cluster data in memory and write
if len(data) > 0:
    comb = combine_cluster(data)
    if args.verbose:
        print("Writing {:d} clusters to run name {:s}"
              .format(comb.id.size, args.outname))
    write_cluster(comb, args.outname, fmt=fmt)

# Step 5: read all integrated data into memory
data = []
for b in basenames:
    in_files = glob.glob(b + "_integrated_*")
    if len(in_files) > 0:
        rinfo={}
        try:
            if args.verbose:
                print("Reading integrated data from run "+b)
            indata = read_integrated(b, nofilterdata=True, read_info=rinfo)
            data.append(indata)
            for k in rinfo.keys():
                if k != 'format':
                    files_read.append(rinfo[k])
                elif fmt is None:
                    fmt = rinfo['format']
        except IOError as e:
            if args.halt:
                raise IOError(e)
            else:
                warnings.warn("IOError during processing of integrated "
                              "data for run " + b +":" 
                              + str(e) + "; processing continues", 
                              RuntimeWarning)

# Step 6: combine integrated data in memory and write
if len(data) > 0:
    comb = combine_integrated(data)
    if args.verbose:
        print("Writing {:d} times, {:d}} trials to run name {:s}"
              .format(comb.time.size, data.target_mass.shape[-1], 
                      args.outname))
    write_integrated(comb, args.outname, fmt=fmt)

# Step 7: read and combine summary files
ntrials = 0
summary_lines = []
read_file = False
for b in basenames:
    try:
        fp = open(b+'_summary.txt', 'r')
        files_read.append(b+'_summary.txt')
        for line in fp:
            if not read_file:
                summary_lines.append(line)
            linesplit = line.split()
            if linesplit[0] == 'n_trials':
                ntrials += int(linesplit[1])
        fp.close()
        read_file = True
    except IOError as e:
        if args.halt:
            raise IOError(e)
        else:
            warnings.warn("IOError during processing of summary file "
                          "for run " + b +":" 
                          + str(e) + "; processing continues", 
                          RuntimeWarning)
if ntrials > 0:
    if args.verbose:
        print("Writing new summary file")
    fpout = open(args.outname+'_summary.txt', 'w')
    for line in summary_lines:
        linesplit = line.split()
        if linesplit[0] == 'model_name':
            fpout.write("model_name           "+osp.basename(args.outname)+"\n")
        elif linesplit[0] == 'out_dir':
            fpout.write("out_dir              "+osp.dirname(args.outname)+"\n")
        elif linesplit[0] == 'n_trials':
            fpout.write("n_trials             {:d}\n".format(ntrials))
        else:
            fpout.write(line)
    fpout.close()

# Step 8: delete files if requested
if args.delete:
    if args.verbose and len(files_read) > 0:
        print("Deleting files:")
    for r in files_read:
        if args.verbose:
            print(r)
        os.remove(r)
