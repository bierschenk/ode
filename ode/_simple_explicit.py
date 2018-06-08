# -*- coding: utf-8 -*-
# Simple explicit integrators
# Explicit integrators calculate the future state from the current state.

__all__ = [
        'verlet', 'iverlet',
          ]


from . _functions import _t_gen


def iverlet(*, ddot_func, x_zero, v_zero, t_range, t_step):
    '''Verlet method integration
    A second order DE of type ddot(x(t)) = A(x(t)) with initial conditions
    x(t_0) = x_0 and dot(x(t_0)) = v_0:
    x1 = x0 + v0*dt + (1/2)*A(x0)*dt^2
    x_{n+1} = 2*x_n - x_{n-1} + A(x_n)*dt^2
    Input: ddot_func, x_zero, v_zero, t_range, t_step
    Output: yields (t, x) until the end of the time range.'''
    t_gen = _t_gen(t_range=t_range, t_step=t_step)
    x = x_zero
    first_step = True
    while True:
        yield (next(t_gen), x)
        if first_step is True:
            x, xold = [sum(i) for i in zip(
                        x,
                        [j * t_step for j in v_zero],
                        [k * t_step * t_step*.5 for k in ddot_func(X=x)]
                        )], x
            first_step = False
        else:
            x, xold = [sum(i) for i in zip(
                        [2 * j for j in x],
                        [-1 * k for k in xold],
                        [d * t_step * t_step for d in ddot_func(X=x)]
                        )], x


def verlet(*, ddot_func, x_zero, v_zero, t_range, t_step):
    '''Verlet method integration
    A second order DE of type ddot(x(t)) = A(x(t)) with initial conditions
    x(t_0) = x_0 and dot(x(t_0)) = v_0:
    x1 = x0 + v0*dt + (1/2)*A(x0)*dt^2
    x_{n+1} = 2*x_n - x_{n-1} + A(x_n)*dt^2
    Input: ddot_func, x_zero, v_zero, t_range, t_step
    Output: (t, x) '''
    t, x = zip(*list(iverlet(
        ddot_func=ddot_func, x_zero=x_zero, v_zero=v_zero,
        t_range=t_range, t_step=t_step)))
    return t, x


def velocity_verlet():
    '''/wiki/Verlet_integration#velocity_verlet
    not quite the same as leapfrog method'''
    pass


def leapfrog():
    '''similar to velocity_verlet but position and velocity are calculated
    at staggered times.'''
    pass
