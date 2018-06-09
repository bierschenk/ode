# -*- coding: utf-8 -*-
# integrators

__all__ = [
        'Euler', 'euler',
        'Verlet', 'verlet',
        'BackwardEuler', 'backwardeuler',
          ]


import math


class Integrator:
    def __iter__(self):
        return self


class ConstantTimestep(Integrator):
    '''The __init__ function of this class sets instance variables for
    integrators with a constant timestep.'''
    def __init__(self, dfun, xzero, timerange, timestep):
        self.timestart, self.timeend = timerange
        assert ((self.timeend - self.timestart) / timestep) > 0, (
            "'timerange' and 'timestep' not consistant. "
            "Check signs and order.")
        self.dfun = dfun
        self.x = xzero
        self.timestep = timestep
        self.time = self.timestart
        self.direction = timestep/abs(timestep)
        self.steps = math.ceil((self.timeend - self.timestart) / timestep)
        self.stepcounter = 0
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
        if self.stepcounter < self.steps:
            if self.status == 'initialized':
                self.status = 'running'
                return self.time, self.x
            else:
                self.stepcounter += 1
                dx = self.dfun(self.time, self.x)
                dxdt = [x * self.timestep for x in dx]
                self.time, self.x = (
                    self.timestart + (self.stepcounter * self.timestep),
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
        super().__init__(dfun, xzero, timerange, timestep)

    def __next__(self):
        if self.stepcounter < self.steps:
            if self.status == 'initialized':
                ddx = self.dfun(self.time, self.x)
                self.xnext = [sum(i) for i in zip(
                        self.x,
                        [v * self.timestep for v in self.v],
                        [(1 / 2) * a * (self.timestep**2) for a in ddx])]
                self.status = 'running'
                return self.time, self.x, self.v
            else:
                self.stepcounter += 1
                ddx = self.dfun(self.time + self.timestep, self.xnext)
                self.time, self.xold, self.x, self.xnext = [
                    self.timestart + (self.stepcounter * self.timestep),
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


class BackwardEuler(ConstantTimestep):
    '''Backward Euler method integration. This class implements a generator.
    https://en.wikipedia.org/wiki/Backward_Euler_method
    Inputs:
        dfun:
            Derivative function of the system.
            The differential system arranged as a series of first-order
            equations: \dot{X} = dfun(t, x)
        xzero:
            The initial condition of the system.
        vzero:
            The initial condition of first derivative of the system.
        timerange:
            The start and end times as (starttime, endtime).
        timestep:
           The timestep.
        convergencethreshold:
            Each step requires an iterative solution of an implicit equation.
            This is the threshold of convergence.
        maxiterations:
            Maximum iterations of the implicit equation before raising
            an exception.
    Outputs:
        t, x:
            for each iteration.
    '''
    def __init__(self, dfun, xzero, timerange, timestep,
                 convergencethreshold=0.0000000001, maxiterations=1000):
        self.convergencethreshold = convergencethreshold
        assert convergencethreshold > 0, 'convergencethreshold must be > 0'
        self.maxiterations = maxiterations
        assert maxiterations > 0, 'maxiterations must be > 0'
        super().__init__(dfun, xzero, timerange, timestep)

    def __next__(self):
        if self.stepcounter < self.steps:
            if self.status == 'initialized':
                self.status = 'running'
                return self.time, self.x
            else:
                self.stepcounter += 1
                iterations = 0
                error = 1 + self.convergencethreshold
                xn1 = self.x
                while (
                        (error >= self.convergencethreshold) and
                        (iterations < self.maxiterations)):
                    iterations += 1
                    xn2 = [i + (j * self.timestep) for i, j in zip(
                        self.x, self.dfun(self.time, xn1))]
                    error = sum([abs(i - j) for i, j in zip(xn1, xn2)])
                    xn1 = xn2
                if error <= self.convergencethreshold:
                    self.time, self.x = (
                        self.timestart + (self.stepcounter * self.timestep),
                        xn1)
                    return self.time, self.x
                else:
                    raise RuntimeError('maximum iterations exceeded')
        else:
            self.status = 'finished'
            raise StopIteration


def backwardeuler(dfun, xzero, timerange, timestep):
    '''Backward Euler method integration.
    This function wraps the Backward Euler class.
    https://en.wikipedia.org/wiki/Backward_Euler_method
    Inputs:
        dfun:
            Derivative function of the system.
            The differential system arranged as a series of first-order
            equations: \dot{X} = dfun(t, x)
        xzero:
            The initial condition of the system.
        vzero:
            The initial condition of first derivative of the system.
        timerange:
            The start and end times as (starttime, endtime).
        timestep:
           The timestep.
        convergencethreshold:
            Each step requires an iterative solution of an implicit equation.
            This is the threshold of convergence.
        maxiterations:
            Maximum iterations of the implicit equation before raising
            an exception.
    Outputs:
        t, x:
    '''
    return zip(*list(BackwardEuler(dfun, xzero, timerange, timestep)))
