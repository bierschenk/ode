# -*- coding: utf-8 -*-
# Simple explicit integrators
# Explicit integrators calculate the future state from the current state.

import itertools
from . _functions import _t_gen



def euler(*, dot_func, x_zero, t_range, t_step):
    '''Euler method integration'''
    def x_n_gen(*, dot_func, x, t_step):
        while True:
            yield x
            dx = [i * t_step for i in dot_func(X = x)] # dx = dt*dot(x)
            x = [sum(j) for j in zip(x, dx)] # x = x + dx
    
    x_n = x_n_gen(dot_func = dot_func,
            x = x_zero, t_step = t_step)
    t = list(_t_gen(
            t_start = t_range[0],
            t_end = t_range[1],
            t_step = t_step,
            ))
    x = list(itertools.islice(x_n, len(t)))
    return t,x


def verlet(*, dot_func, x, t_step):
    '''Verlet method integration
    A second order DE of type ddot(x(t)) = A(x(t)) with initial conditions
    x(t_0) = x_0 and dot(x(t_0)) = v_0:
    x1 = x0 + v0*dt + (1/2)*A(x0)*dt^2
    x_{n+1} = 2*x_n - x_{n-1} + A(x_n)*dt^2'''
    pass


def velocity_verlet():
    '''/wiki/Verlet_integration#velocity_verlet
    not quite the same as leapfrog method'''
    pass


def leapfrog():
    '''similar to velocity_verlet but position and velocity are calculated
    at staggered times.'''
    pass



#change euler to two functions: ieuler and euler. ieuler returns a generator, euler 
# calls list(ieuler).

