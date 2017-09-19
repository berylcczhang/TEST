"""
This defines a class that can be used to draw random numbers from a
normal distribution.
"""

import numpy as np
import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from scipy.special import erf
else:
    # Dummy erf function
    def erf(x):
        pass
from .slug_pdf_segment import slug_pdf_segment

class slug_pdf_normal(slug_pdf_segment):
    """
    This class defines a normal segment of a PDF
    """

    def __init__(self, a, b, rand=None, fp=None, mean=None, disp=None):
        """
        Class initializer for a PDF of the form 
        dp/dx ~ e^(-(x - mean)^2 / (2 sigma^2) in the range
        [a,b]

        Parameters
           a : float
              lower limit of the segment
           b : float
              upper limit of the segment
           rand : RandomState
              a numpy RandomState object, used to produce random
              deviates; if None, a new RandomState object is created
           mean : float
              location of the peak of the PDF; ignored if fp is set
           disp : float
              dispersion of the PDFl; ignored if fp is set
           fp : file
              file object pointing to the start of the PDF data
        """

        # Call parent constructor
        super(slug_pdf_normal, self).__init__(a, b, rand)

        # See if we were given a file or explicit values
        if fp is None:
            if mean is None or disp is None:
                raise ValueError( 
                    "slug_pdf_normal.__init__: "+
                    "must set either fp or both mean and disp")
            else:
                self._mean = mean
                self.disp = disp

        else:
            # Read mean and dispersion from the file we've been given
            self._mean = None
            self._disp = None
            for line in fp:

                # Strip trailing comments
                while line.rfind('#') != -1:
                    line = line[:line.rfind('#')]

                # Strip leading and trailing whitespace
                line = line.strip()

                # Skip blank lines
                if len(line) == 0:
                    continue

                # Read the keyword and store its value
                spl = line.split()
                if spl[0] == 'mean':
                    if len(spl) != 2:
                        raise IOError(errno.EIO, 
                            "slug_pdf_lognormal: expected "+
                            "'mean MEAN', found "+line)
                    self._mean = float(spl[1])
                elif spl[0] == 'disp':
                    if len(spl) != 2:
                        raise IOError(errno.EIO, 
                            "slug_pdf_lognormal: expected "+
                            "'disp DISP', found "+line)
                    self._disp = float(spl[1])

                # Break if we're done
                if self._mean is not None and self._disp is not None:
                    break

            # Normalize
            self.normalize()


    @slug_pdf_segment.a.setter
    def a(self, val):
        self._a = val
        self.normalize()
    @slug_pdf_segment.b.setter
    def b(self, val):
        self._b = val
        self.normalize()
    @property
    def mean(self):
        return self._mean
    @mean.setter
    def mean(self, val):
        self._mean = val
        self.normalize()
    @property
    def disp(self):
        return self._disp
    @disp.setter
    def disp(self, val):
        self._disp = val
        self.normalize()

        
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
            pdf = self._norm * \
                  np.exp(-(xarr-self._mean)**2 /
                           (2.0*self._disp**2))
            pdf[np.logical_or(xarr < self._a, xarr > self._b)] = 0.0
        else:
            if self._a <= x and x <= self._b:
                pdf = self._norm * \
                      np.exp(-(x-self._mean)**2 /
                               (2.0*self._disp**2))
            else:
                pdf = 0.0
        return pdf
        

    def draw(self, *d):
        """
        Draw from the lognormal PDF

        Parameters:
           d0, d1, ..., dn : int, optional
              Dimensions of the returned array; if left unspecified, a
              single python float is returned

        Returns:
           x : float or array
              One or more numbers drawn from the PDF
        """

        # Scalar case
        if len(d) == 0:
            while True:
                sample = self._rnd.normal(
                    mean=np.log(self._mean), scale=self._disp)
                if self._a <= sample and sample <= self._b:
                    return sample

        # Vector case
        samples = np.zeros(*d)
        samples[:] = np.nan
        while True:
            idx = np.isnan(samples)
            nidx = np.sum(idx)
            if nidx == 0:
                return samples
            s = self._rng.normal(
                mean=np.log(self._mean), sigma=self._disp,
                size=nidx)
            accept = np.logical_and(s >= self._a, s <= self._b)
            samples[np.where(idx)[0][accept]] = s[accept]


    def expectation(self):
        """
        Compute the expectation value of the lognormal PDF

        Parameters:
           None

        Returns:
           expec : float
              Expectation value of the PDF
        """
        return 2.0*self._disp*(self._a-self._b)/self._norm + \
            np.sqrt(2.0*np.pi)*self._mean * \
            (erf((self._b-self._mean)/(np.sqrt(2.0)*self._disp)) -
             erf((self._a-self._mean)/(np.sqrt(2.0)*self._disp)))
    

    def normalize(self):
        """
        Recompute the normalization factor for the PDF

        Parameters:
           None

        Returns:
           Nothing
        """
        self._norm \
            = np.sqrt(2.0/np.pi) / self._disp / \
            (erf((self._b-self._mean)/(np.sqrt(2.0)*self._disp)) -
             erf((self._a-self._mean)/(np.sqrt(2.0)*self._disp)))
