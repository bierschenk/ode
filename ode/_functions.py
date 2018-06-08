# -*- coding: utf-8 -*-
# helper functions

__all__ = [
        '_t_gen',
          ]
import math


def _t_gen(t_range, t_step):
    t_start, t_end = t_range
    t = t_start
    n = math.ceil((t_end - t_start) / t_step)
    for i in range(1, n+2):
        yield t
        t = t_start + i*t_step


def _newton_raphson():
    pass
