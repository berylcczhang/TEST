"""
This defines a class that can be used to draw random numbers from a
delta distribution.
"""

import numpy as np
from .slug_pdf_segment import slug_pdf_segment

class slug_pdf_delta(slug_pdf_segment):
    """
    This class defines a delta function segment of a PDF
    """

    def __init__(self, a):
        """
        Class initializer for a PDF of the form dp/dx ~ delta(x-a)

        Parameters
            a : float
        """

        # Call parent constructor, setting rand to 0 instead of None so
        # that we don't generate an unnecessary RNG
        super(slug_pdf_delta, self).__init__(a, a, rand=0)
      

    # Setters; note that delta PDFs have only a single value, so a and
    # b are always the same
    @slug_pdf_segment.a.setter
    def a(self, val):
      self._a = val
      self._b = val
    @slug_pdf_segment.b.setter
    def b(self, val):
      self._a = val
      self._b = val

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
            pdf[xarr == self._a] == np.inf
        else:
            if x == self._a:
              pdf = np.inf
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
        if len(d) == 0:
          return self._a
        else:
          x = np.zeros(*d)
          x[:] = self._a
          return x


    def expectation(self):
        """
        Compute the expectation value of the lognormal PDF

        Parameters:
           None

        Returns:
           expec : float
              Expectation value of the PDF
        """
        return self._a
