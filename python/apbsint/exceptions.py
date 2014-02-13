"""
exceptions
==========

Exception classes for ApBsInT modules.

"""

__all__ = ['ApBsWrapError']

# Exception classes

class ApBsWrapError(Exception):
    """
    ApBsWrapError
    =============

    Exception raised in C++ code.

    Attributes:
        msg -- explanation of the error

    """
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg
