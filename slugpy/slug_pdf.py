"""
This defines a class that can be used to parse slug PDF files and draw
random values from them. This code is thread-safe, in the sense that
if multiple slug_pdf instances are instantiated in different threads,
the random streams they generate will not be identical.
"""

# Special check for readthedocs
import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from numpy.random import RandomState
else:
    def RandomState():
        return None

import numpy as np
import errno
from .slug_pdf_delta import slug_pdf_delta
from .slug_pdf_exponential import slug_pdf_exponential
from .slug_pdf_powerlaw import slug_pdf_powerlaw
from .slug_pdf_lognormal import slug_pdf_lognormal
from .slug_pdf_powerlaw import slug_pdf_powerlaw
from .slug_pdf_normal import slug_pdf_normal

class slug_pdf(object):
    """
    A class that implements the SLUG PDF drawing method. This class
    contains a method to parse slug-formatted PDF files and then draw
    values from the PDFs they specify. This class is thread-safe, in
    the sense that if multiple slug_pdf instances are instantiated in
    different threads, the random streams they generate will not be
    identical.
    """

    ##################################################################
    # Initialization method
    ##################################################################
    def __init__(self, pdffile=None):
        self.segments = []
        self.wgts = np.array(0)
        # Use /dev/urandom to initialize, so should be thread-safe
        self._rnd = RandomState()
        if pdffile is not None:
            self.readfile(pdffile)

    ##################################################################
    # Method to parse PDF files
    ##################################################################
    def readfile(self, pdffile):

        # Open the file
        fp = open(pdffile, 'r')

        # Read until we see breakpoints or advanced
        for line in fp:

            # Strip trailing comments
            while line.rfind('#') != -1:
                line = line[:line.rfind('#')]

            # Strip leading and trailing whitespace
            line = line.strip()

            # Skip blank lines
            if len(line) == 0:
                continue

            # Is this line advanced or breakpoints? If not, raise
            # error.
            if line.split()[0] == 'breakpoints':
                basic = True
            elif line.strip() == 'advanced':
                basic = False
            else:
                raise IOError(errno.EIO, "slug_pdf: expected 'breakpoints' "+
                              "or 'advanced', found '"+line+"'")

            # If we are in basic mode, extract the breakpoints
            if basic:
                if len(line.split()) < 3:
                    raise IOError(errno.EIO, 
                        "slug_pdf: need at least 2 breakpoints")
                self.bkpts = np.array(line.split()[1:], dtype='d')

            # Done reading header
            break

        # Make sure we haven't reached EOF with setting a mode
        if basic is None:
            raise IOError(errno.EIO, "slug_pdf: didn't find 'breakpoints' "+
                          "or 'advanced' in file "+pdffile)
        
        # Now parse the rest of the file
        if basic:
            self.parse_basic(fp)
        else:
            self.parse_advanced(fp)

        # Close file
        fp.close()

    ##################################################################
    # Method to parse basic PDF files
    ##################################################################
    def parse_basic(self, fp):

        # Flag if we're in a segment
        in_segment = False

        # Initialize segment list
        self.segments = []
        
        # Loop through file
        for line in fp:

            # Strip trailing comments
            while line.rfind('#') != -1:
                line = line[:line.rfind('#')]

            # Strip leading and trailing whitespace
            line = line.strip()

            # Skip blank lines
            if len(line) == 0:
                continue

            # Parse command
            if line == 'segment':

                # We're being asked to start a new segment
                if not in_segment:
                    in_segment = True
                else:
                    raise IOError(errno.EIO, "slug_pdf: expected 'type "+
                                  "TYPENAME', found"+
                                  " '"+line+"'")

            elif line.split()[0] == 'type':

                # We're being given a segment type; make sure we're in
                # a segment
                if not in_segment:
                    raise IOError(errno.EIO, "slug_pdf: expected 'segment',"+
                                  " found '"+line+"'")

                # Make sure we have a valid segment type
                if len(line.split()) != 2:
                    raise IOError(errno.EIO, "slug_pdf: expected 'type "+
                                  "TYPENAME', found"+
                                  " '"+line+"'")

                # Read segment type
                segtype = line.split()[1]

                # Grab breakpoints for this segment
                ptr = len(self.segments)
                if ptr > len(self.bkpts)-1:
                    raise IOError(errno.EIO, 
                        "slug_pdf: expected {:d}".len(self.bkpts-1)+
                        " segments, found too many")
                a = self.bkpts[ptr]
                b = self.bkpts[ptr+1]
                
                # Call parsing routine for this segment
                if segtype == 'delta':
                    if a != b:
                        raise IOError(errno.EIO, "slug_pdf: delta segments "+
                                      "must have a == b")
                    self.segments.append(slug_pdf_delta(a))
                elif segtype == 'exponential':
                    self.segments.append(
                        slug_pdf_exponential(a, b, fp=fp, rand=self._rnd))
                elif segtype == 'lognormal':
                    self.segments.append(
                        slug_pdf_lognormal(a, b, fp=fp, rand=self._rnd))
                elif segtype == 'powerlaw':
                    self.segments.append(
                        slug_pdf_powerlaw(a, b, fp=fp, rand=self._rnd))
                elif segtype == 'normal':
                    self.segments.append(
                        slug_pdf_normal(a, b, fp=fp, rand=self._rnd))
                elif segtype == 'schechter':
                    raise NotImplementedError(
                        "slug_pdf: schechter segments not yet "+
                        "implemented")

                # Reset the segment pointer
                in_segment = False

        # Make sure we got enough segments
        if len(self.segments) != len(self.bkpts)-1:
            raise IOError(errno.EIO, 
                "slug_pdf: expected {:d} segments, found {:d}".
                format(len(self.bkpts)-1, len(self.segments)))

        # Convert segments to array
        self.segments = np.array(self.segments)
        
        # Compute the relative weights required to make the segments
        # continuous
        self.wgts = np.zeros(len(self.segments))
        self.wgts[0] = 1.0
        for i in range(1, len(self.segments)):
            self.wgts[i] = self.wgts[i-1] * \
                           self.segments[i-1](self.segments[i-1].b) / \
                           self.segments[i](self.segments[i].a)

        # Set weights to give overall normalization of 1
        self.wgts = self.wgts / np.sum(self.wgts)

        # Store cumulative weights
        self.cumwgts = np.cumsum(self.wgts)

        
    ##################################################################
    # Method to parse advanced PDF files
    ##################################################################
    def parse_advance(self, fp):
        raise NotImplementedError( 
            "slug_pdf: advanced mode files not yet supported")

    ##################################################################
    # Methods to evaluate PDF and to draw from it
    ##################################################################
    def __call__(self, x):
        """
        Return the value of the PDF evaluated at x

        Parameters:
           x : float or array
              The value(s) at which to evaluate the PDF

        Returns:
           pdf : float or array
              The value of the PDF evaluated at x
        """
        if hasattr(x, '__iter__'):
            xarr = np.array(x)
            pdf = np.zeros(xarr.shape)
            for s, w in zip(self.segments, self.wgts):
                pdf += w*s(xarr)
        else:
            pdf = 0.0
            for s, w in zip(self.segments, self.wgts):
                pdf += w*s(x)
        return pdf

    def draw(self, *d):
        """
        Draw from the PDF

        Parameters:
           d0, d1, ..., dn : int, optional
              Dimensions of the returned array; if left unspecified, a
              single python float is returned

        Returns:
           x : float or array
              One or more numbers drawn from the PDF
        """

        # Figure out which segments to draw from
        dev = self._rnd.rand(*d)
        seg = np.digitize(dev, self.cumwgts)

        # Draw
        if hasattr(dev, '__iter__'):
            samples = np.zeros(dev.shape)
            for i, s in enumerate(self.segments):
                idx = (i == seg)
                nidx = np.sum(idx)
                if nidx > 0:
                    samples[idx] = s.draw(nidx)
        else:
            samples = self.segments[seg].draw()

        # Return
        return samples
