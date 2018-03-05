# -*- coding: utf-8 -*-
# Simple implicit integrators
# Implicit integrators calculate the future state from the future state
# according to y_{k+1} = y_k + h*f(t_{k+1}, y_{k_1}) where f(t,y) = dy/dt.


import itertools
from . _functions import _t_gen


def ibackward_euler(*, dot_func, x_zero, t_range, t_step,
        convergence_threshold = 0.0000000001, max_iterations = 1000):
    '''Backward Euler, or Implicit Euler method integration as a generator.
    X_{n+1} = X_n + dt*dot(X_{n+1})
    Iteration is used to find X_{n+1}. First iteration is Euler
    (X_n + dt*dot(X_n)).
    Input: derivative function, X_zero vector, time range, time step
    Output: yeilds (t,x) until end of time range. x is the state vector
    as a tuple.'''
    t_gen = _t_gen(t_range = t_range, t_step = t_step)
    x = x_zero
    while True:
        yield (next(t_gen), x)
        x1_i = x
        error = convergence_threshold+1
        converged = False
        iterations = 0
        while not converged:
            dx = [i * t_step for i in dot_func(X = x1_i)] # dx = dt*dot(x1_i)
            x1_i2 = [sum(j) for j in zip(x,dx)] # this is x1_i2 = x + dx
            error = sum([abs(i[0] - i[1]) for i in zip(x1_i, x1_i2)])
            iterations +=1
            if error <= convergence_threshold:
                x = x1_i2
                converged = True
            elif iterations > max_iterations:
                x = 'Convergence Error'
                converged = True
            else:
                x1_i = x1_i2


def backward_euler(*, dot_func, x_zero, t_range, t_step,
        convergence_threshold = 0.0000000001, max_iterations = 1000):
    '''Backward Euler, or Implicit Euler method integration.
    X_{n+1} = X_n + dt*dot(X_{n+1})
    Iteration is used to find X_{n+1}. First iteration is Euler
    (X_n + dt*dot(X_n)).'''
    t, x = zip(*list(ibackward_euler(dot_func = dot_func, x_zero = x_zero,
            t_range = t_range, t_step = t_step,
            convergence_threshold = 0.0000000001, max_iterations = 1000)))
    return t, x
