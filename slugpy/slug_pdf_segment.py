"""
This defines a class of a single segment of a PDF; used together with slug_pdf.
"""

# Special check for readthedocs
import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from numpy.random import RandomState
else:
    def RandomState():
        return None

class slug_pdf_segment(object):
    """
    A class that works together with slug_pdf to implement the PDF
    drawing method used by slug. This is a purely abstract class, used
    to define a common interface for PDF segments.
    """
    
    def __init__(self, a, b, rand=None):
        """
        Initializer for slug_pdf_segment

        Parameters
           a : float
              lower limit of the segment
           b : float
              upper limit of the segment
           rand : RandomState
              a numpy RandomState object, used to produce random
              deviates; handled in user-space for thread safety
              reasons
        """

        # Store limits and RNG
        self._a = a
        self._b = b
        if rand is not None:
            self._rnd = rand
        else:
            self._rnd = RandomState()

    def draw(self, *d):
        raise NotImplementedError(1, 
            "slug_pdf_segment.draw() should never be invoked directly!")

    def expectation(self):
        raise NotImplementedError(1, 
            "slug_pdf_segment.expectation() should never be invoked"
            "directly!")

    def __call__(self):
        raise NotImplementedError(1, 
            "slug_pdf_segment() should never be invoked directly!")

    @property
    def a(self):
        return self._a
    @a.setter
    def a(self, val):
        self._a = val
    @property
    def b(self):
        return self._b
    @b.setter
    def b(self, val):
        self._b = val

