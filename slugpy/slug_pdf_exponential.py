"""
This defines a class that can be used to draw random numbers from an
exponential distribution.
"""

import numpy as np
from .slug_pdf_segment import slug_pdf_segment

class slug_pdf_exponential(slug_pdf_segment):
    """
    This class defines an exponential segment of a PDF
    """

    def __init__(self, a, b, rand=None, fp=None, scale=None):
        """
        Class initializer for a PDF of the form 
        dp/dx ~ e^(x/scale) in the range [a,b]

        Parameters
           a : float
              lower limit of the segment
           b : float
              upper limit of the segment
           rand : RandomState
              a numpy RandomState object, used to produce random
              deviates; if None, a new RandomState object is created
           scale : float
              scale length of the exponential; ignored if fp is set
           fp : file
              file object pointing to the start of the PDF data
        """

        # Call parent constructor
        super(slug_pdf_exponential, self).__init__(a, b, rand)

        # See if we were given a file or explicit values
        if fp is None:
            if scale is None:
                raise ValueError( 
                    "slug_pdf_exponential.__init__: "+
                    "must set either fp or scale")
            else:
                self.scale = scale

        else:
            # Read mean and dispersion from the file we've been given
            self._scale = None
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
                if spl[0] == 'scale':
                    if len(spl) != 2:
                        raise IOError(errno.EIO, 
                            "slug_pdf_exponential: expected "+
                            "'scale SCALE', found "+line)
                    self.scale = float(spl[1])
                    break


    @slug_pdf_segment.a.setter
    def a(self, val):
        self._a = val
        self.normalize()
    @slug_pdf_segment.b.setter
    def b(self, val):
        self._b = val
        self.normalize()
    @property
    def scale(self):
        return self._scale
    @scale.setter
    def scale(self, val):
        self._scale = val
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
            pdf = self._norm * np.exp(-x/self._xscale)
            pdf[np.logical_or(xarr < self._a, xarr > self._b)] = 0.0
        else:
            if self._a <= x and x <= self._b:
                pdf = self._norm * np.exp(x/self._xscale)
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

        dev = self._rnd.rand(*d)
        samples = -self._scale * \
                  ( np.log(dev*np.exp(-self._b/self._scale) +
                           (1.0-dev)*np.exp(-self._a/self._scale)) )


    def expectation(self):
        """
        Compute the expectation value of the lognormal PDF

        Parameters:
           None

        Returns:
           expec : float
              Expectation value of the PDF
        """
        return self._a + self._scale + (self._a-self._b) * \
            np.exp(self._a/self._scale) / \
            (np.exp(self._b/self._scale) - np.exp(self._a/self._scale))
    

    def normalize(self):
        """
        Recompute the normalization factor for the PDF

        Parameters:
           None

        Returns:
           Nothing
        """
        self._norm \
            = 1.0 / (self._scale *
                     (np.exp(-self._a/self._scale) -
                      np.exp(-self._b/self._scale)))
