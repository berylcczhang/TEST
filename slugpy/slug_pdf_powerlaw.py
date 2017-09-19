"""
This defines a class that can be used to draw random numbers from a
powerlaw distribution.
"""

import numpy as np
from .slug_pdf_segment import slug_pdf_segment

class slug_pdf_powerlaw(slug_pdf_segment):
    """
    This class defines a powerlaw segment of a PDF
    """

    def __init__(self, a, b, rand=None, fp=None, alpha=None):
        """
        Class initializer for a PDF of the form x^alpha in the range
        [a,b]

        Parameters
           a : float
              lower limit of the segment
           b : float
              upper limit of the segment
           rand : RandomState
              a numpy RandomState object, used to produce random
              deviates; if None, a new RandomState object is created
           alpha : float
              slope of the PDF; ignored if fp is set
           fp : file
              file object pointing to the start of the PDF data
        """

        # Call parent constructor
        super(slug_pdf_powerlaw, self).__init__(a, b, rand)

        # See if we were given a file or an explicit slope
        if fp is None:
            if alpha is None:
                raise ValueError( 
                    "slug_pdf_powerlaw.__init__: "+
                    "must set fp or alpha")
            else:
                self.alpha = alpha

        else:
            # Read slope from the file we've been given
            for line in fp:

                # Strip trailing comments
                while line.rfind('#') != -1:
                    line = line[:line.rfind('#')]

                # Strip leading and trailing whitespace
                line = line.strip()

                # Skip blank lines
                if len(line) == 0:
                    continue

                # Make sure this is a properly formatted slope specified
                if len(line.split()) != 2 or line.split()[0] != 'slope':
                    raise IOError(errno.EIO, 
                        "slug_pdf_powerlaw: expected 'slope SLOPE', "+
                        "found '"+line+"'")

                # Read slope
                self.alpha = float(line.split()[1])

                # Break
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
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, val):
        self._alpha = val
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
            pdf = self._norm * xarr**self._alpha
            pdf[np.logical_or(xarr < self._a, xarr > self._b)] = 0.0
        else:
            if self._a <= x and x <= self._b:
                pdf = self._norm * x**self._alpha
            else:
                pdf = 0.0
        return pdf
        

    def draw(self, *d):
        """
        Draw from the powerlaw PDF

        Parameters:
           d0, d1, ..., dn : int, optional
              Dimensions of the returned array; if left unspecified, a
              single python float is returned

        Returns:
           x : float or array
              One or more numbers drawn from the PDF
        """

        dev = self._rnd.rand(*d)
        if self._alpha != -1.0:
            return (dev*self._b**(self._alpha+1) +
                    (1-dev)*self._a**(self._alpha+1)) \
                    **(1.0/(self._alpha+1.0))
        else:
            return self._b**dev / self._a**(dev-1.0)


    def expectation(self):
        """
        Compute the expectation value of the lognormal PDF

        Parameters:
           None

        Returns:
           expec : float
              Expectation value of the PDF
        """
        if (self._alpha != -1.0) and (self._alpha != -2.0):
            return (self._alpha+1.0)/(self._alpha+2.0) * \
                (self._b**(self._alpha+2.0) -
                 self._a**(self._alpha+2.0)) / \
                 (self._b**(self._alpha+1.0) -
                  self._a**(self._alpha+1.0))
        elif self._alpha == -1.0:
            return (self._b-self._a) / np.log(self._b/self._a)
        else:
            return (self._b+self._a)/2.0

    def normalize(self):
        """
        Recompute the normalization factor for the PDF

        Parameters:
           None

        Returns:
           Nothing
        """
        if self._alpha != -1.0:
            self._norm = (self.alpha+1.0) / \
                         (self._b**(self._alpha+1.0) - 
                          self._a**(self._alpha+1.0))
        else:
            self._norm = 1.0 / np.log(self._b/self._a)
            
            
