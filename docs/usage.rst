=====
Usage
=====

To use ode in a project::

    import ode

The following integration methods are included in ode:
 * Euler's method
 * Backward Euler method

The integration methods use a derivative function (*dot_func*) to integrate.
The system is factored into a series of first-order differential equations.
The derivative function is derived from the system parameters. It returns the 
derivative of the system calculated from the current state of the system.

The current state and derivative of the system are represented as lists.

Euler's method
-----

X_{n+1} = X_n + dt*dot(X_n)

Backward Euler method
-----

X_{n+1} = X_n + dt*dot(X_{n+1})
