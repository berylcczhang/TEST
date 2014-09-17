"""
Routine to read in a SLUG summary output file.
"""

import os
import os.path as osp

def read_summary(model_name, output_dir=None):
    """
    Function to open a SLUG output summary file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set

    Returns
       summary : dict
          A dict containing all the keywords stored in the output file

    Raises
       IOError, if a summary file for the specified model cannot be found
    """

    # Did we get a specific directory in which to look? If not, try
    # current directory
    if output_dir is None:
        outdir = "."
    else:
        outdir = output_dir

    # Try to open in cwd
    fname = osp.join(outdir, model_name+'_summary.txt')
    try:
        fp = open(fname, 'r')
    except IOError:
        fp = None

    # If that failed, and we didn't get an explicit directory
    # specification, try looking in SLUG_DIR/output
    if (fp is None) and (output_dir is None) and \
       ('SLUG_DIR' in os.environ):
        outdir = osp.join(os.environ['SLUG_DIR'], 'output')
        fname = osp.join(outdir, model_name+'_summary.txt')
        fp = open(fname, 'r')

    # Burn first line
    fp.readline()

    # Read data
    summary = {}
    for line in fp:
        linesplit = line.split()
        if linesplit[0] != 'phot_bands':
            # Try converting to integer or float; if that fails, store
            # as string
            try:
                summary[linesplit[0]] = int(linesplit[1])
            except ValueError:
                try:
                    summary[linesplit[0]] = float(linesplit[1])
                except ValueError:
                    summary[linesplit[0]] = linesplit[1]
        else:
            phot = linesplit[1:]
            for p in phot:
                p = p[:-1]   # Strip final comma
            summary[linesplit[0]] = phot

    # Close file and return
    fp.close()
    return summary
