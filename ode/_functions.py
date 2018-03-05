# -*- coding: utf-8 -*-
# helper functions



def _t_gen(*, t_range, t_step):
    t_start, t_end = t_range
    t = t_start
    n = 0
    while t < t_end + t_step:
        yield t
        n += 1
        t  = t_start + (n * t_step)


def _newton_raphson():
    pass
