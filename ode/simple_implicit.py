# -*- coding: utf-8 -*-
# Simple implicit integrators
# Implicit integrators calculate the future state from the future state
# according to y_{k+1} = y_k + h*f(t_{k+1}, y_{k_1}) where f(t,y) = dy/dt.


import itertools
from . _functions import _t_gen



def backward_euler(*, dot_func, x_zero, t_range, t_step,
        convergence_threshold = 0.0000000001, max_iterations = 1000):
    '''Backward Euler, or Implicit Euler method integration.
    X_{n+1} = X_n + dt*dot(X_{n+1})
    Iteration is used to find X_{n+1}. First iteration is Euler
    (X_n + dt*dot(X_n)).'''
    def x_n_gen(*, dot_func, x, t_step, convergence_threshold):
        while True:
            yield x
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
    
    x_n = x_n_gen(dot_func = dot_func, x = x_zero, t_step = t_step,
            convergence_threshold = convergence_threshold)
    t = list(_t_gen(
        t_start = t_range[0],
        t_end = t_range[1],
        t_step = t_step,
        ))
    x = list(itertools.islice(x_n, len(t)))
    return t,x
# is it possible to **NOT** return a list. Instead create a generator and yield **NEXT**






