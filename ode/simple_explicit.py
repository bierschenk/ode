# -*- coding: utf-8 -*-
# Simple explicit integrators
# Explicit integrators calculate the future state from the current state.

import itertools
from . _functions import _t_gen


def ieuler(*, dot_func, x_zero, t_range, t_step):
    '''Euler method integration, as a python generator.
    Input: derivative function, X_zero vector, time range, time step
    Output: yeilds (t,x) until end of time range. x is the state vector
    as a tuple.'''
    t_gen = _t_gen(t_range = t_range, t_step = t_step)
    x = x_zero
    while True:
        yield (next(t_gen), x)
        dx = [i * t_step for i in dot_func(X = x)] # dx = dt*dot(x)
        x = [sum(j) for j in zip(x, dx)] # x = x + dx

def euler(*, dot_func, x_zero, t_range, t_step):
    '''Euler method integration'''
    t, x = zip(*list(ieuler(dot_func = dot_func, x_zero = x_zero,
                  t_range = t_range, t_step = t_step)))
    return t, x


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
