# -*- coding: utf-8 -*-
# integrators

__all__ = [
        'Euler', 'euler'
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
