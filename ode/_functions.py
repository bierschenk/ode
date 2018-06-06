# -*- coding: utf-8 -*-
# helper functions

import math
# def _t_gen(t_range, t_step):
#     t_start, t_end = t_range
#     t = t_start
#     n = 0
#     while t < t_end + t_step:
#         yield t
#         n += 1
#         t = t_start + (n * t_step)


def _t_gen(t_range, t_step):
    t_start, t_end = t_range
    t = t_start
    n = math.ceil((t_end - t_start) / t_step)
    for i in range(1, n+2):
        yield t
        t = t_start + i*t_step


def _newton_raphson():
    pass
