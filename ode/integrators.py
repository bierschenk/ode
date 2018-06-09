# -*- coding: utf-8 -*-
# integrators

__all__ = [
        'Euler', 'euler',
        'Verlet', 'verlet',
          ]


# from ._functions import *


class Integrator:
    def __iter__(self):
        return self


class ConstantTimestep(Integrator):
    '''The __init__ function of this class sets instance variables for
    integrators with a constant timestep.'''
    def __init__(self, dfun, xzero, timerange, timestep):
        self.timestart, self.timeend = timerange
        if not ((self.timeend - self.timestart) / timestep) > 0:
            raise ValueError("'timerange' and 'timestep' not consistant."
                             "Check signs and order.")
        self.dfun = dfun
        self.x = xzero
        self.timestep = timestep
        self.time = self.timestart
        self.direction = timestep/abs(timestep)
        self.status = 'initialized'


class Euler(ConstantTimestep):
    '''Euler method integration. This class implements a generator.
    https://en.wikipedia.org/wiki/Euler_method
    Inputs:
        dfun:
            derivative function of the system.
            The differential system arranged as a series of first-order
            equations: \dot{X} = dfun(t, x)
        xzero:
            the initial condition of the system
        timerange:
            the start and end times as (starttime, endtime)
        timestep:
            the timestep
    Outputs:
        t, x:
            for each iteration.
    '''

    def __next__(self):
        if (self.time * self.direction) < self.timeend:
            if self.status == 'initialized':
                self.status = 'running'
                return self.time, self.x
            else:
                dx = self.dfun(self.time, self.x)
                dxdt = [x * self.timestep for x in dx]
                self.time, self.x = (
                    self.time + self.timestep,
                    [xi + didt for xi, didt in zip(self.x, dxdt)])
                return self.time, self.x
        else:
            self.status = 'finished'
            raise StopIteration


def euler(dfun, xzero, timerange, timestep):
    '''Euler method integration. This function wraps the Euler class.
    https://en.wikipedia.org/wiki/Euler_method
    Inputs:
        dfun:
            derivative function of the system.
            The differential system arranged as a series of first-order
            equations: \dot{X} = dfun(t, x)
        xzero:
            the initial condition of the system
        timerange:
            the start and end times as (starttime, endtime)
        timestep:
            the timestep
    Outputs:
        t, x:
    '''
    return zip(*list(Euler(dfun, xzero, timerange, timestep)))


class Verlet(ConstantTimestep):
    '''Verlet method integration. This class implements a generator.
    https://en.wikipedia.org/wiki/Verlet_integration
    Inputs:
        dfun:
            second derivative function of the system.
            The differential system arranged as a series of second-order
            equations: \ddot{X} = dfun(t, x)
        xzero:
            the initial condition of the system
        vzero:
            the initial condition of first derivative of the system
        timerange:
            the start and end times as (starttime, endtime)
        timestep:
            the timestep
    Outputs:
        t, x, v:
            for each iteration.
    '''
    def __init__(self, dfun, xzero, vzero, timerange, timestep):
        self.v = vzero
        self.xold = xzero
        self.xold = xzero
        super().__init__(dfun, xzero, timerange, timestep)

    def __next__(self):
        if (self.time * self.direction) < self.timeend:
            if self.status == 'initialized':
                ddx = self.dfun(self.time, self.x)
                self.xnext = [sum(i) for i in zip(
                        self.x,
                        [v * self.timestep for v in self.v],
                        [(1 / 2) * a * (self.timestep**2) for a in ddx])]
                self.status = 'running'
                return self.time, self.x, self.v
            elif self.status == 'running':
                ddx = self.dfun(self.time + self.timestep, self.xnext)
                self.time, self.xold, self.x, self.xnext = [
                    self.time + self.timestep,
                    self.x,
                    self.xnext,
                    [sum(i) for i in zip(
                        [2 * x for x in self.xnext],
                        [-1 * x for x in self.x],
                        [a * (self.timestep**2) for a in ddx]
                        )]
                    ]
                self.v = [(xnext - xold) / (2 * self.timestep)
                          for xnext, xold in zip(self.xnext, self.xold)]
                return self.time, self.x, self.v
        else:
            self.status = 'finished'
            raise StopIteration


def verlet(dfun, xzero, vzero, timerange, timestep):
    '''Verlet method integration. This function wraps the Verlet class.
    https://en.wikipedia.org/wiki/Verlet_integration
    Inputs:
        dfun:
            second derivative function of the system.
            The differential system arranged as a series of second-order
            equations: \ddot{X} = dfun(t, x)
        xzero:
            the initial condition of the system
        vzero:
            the initial condition of first derivative of the system
        timerange:
            the start and end times as (starttime, endtime)
        timestep:
            the timestep
    Outputs:
        t, x, v:
    '''
    return zip(*list(Verlet(dfun, xzero, vzero, timerange, timestep)))
