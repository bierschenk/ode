# -*- coding: utf-8 -*-

__author__ = """bierschenk"""
__email__ = 'bierschenk.devel@gmail.com'
__version__ = '0.1.0'


# Import all functions transparently to ode.

from . simple_explicit import euler
from . simple_implicit import backward_euler

