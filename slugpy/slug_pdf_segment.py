"""
This defines a class of a single segment of a PDF; used together with slug_pdf.
"""

class slug_pdf_segment(object):
    """
    A class that works together with slug_pdf to implement the PDF
    drawing method used by slug. This is a purely abstract class, used
    to define a common interface for PDF segments.
    """
    
    def __init__(self, a, b, fp=None):
        raise NotImplementedError(
            "slug_pdf_segment should never be instantiated directly!")

    def draw(self, *d):
        raise NotImplementedError(
            "slug_pdf_segment.draw() should never be invoked directly!")

    def expectation(self):
        raise NotImplementedError(
            "slug_pdf_segment.expectation() should never be invoked"
            "directly!")

    def __call__(self):
        raise NotImplementedError(
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

