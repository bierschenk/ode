Integration Methods
===================

The following integration methods are included in ode:
 * Euler's method
 * Backward Euler method
 * Verlet method

The integration methods operate on systems of either first or second
order differential equations. By convention :math:`X` is the vector
containing the state variables of the system, :math:`f(t,X)` is
a function returning either the first or second derivative of the
system, and :math:`h` is the timestep.

The current state and derivative of the system are represented as lists.

Euler's method
--------------
Euler's method is an explicit method for solving a system of
first order differential equations.

.. math::
    f(t,X) = \dot{X}

.. math::
    X_{n+1} = X_n + h \cdot f(t_n, X_n)

Backward Euler method
---------------------
Euler's method is an implicit method for solving a system of
first order differential equations.

.. math::
    f(t,X) = \dot{X}

.. math::
    X_{n+1} = X_n + h \cdot f(t_{n+1}, X_{n+1})

Verlet method
-------------
The Verlet method, also called Störmer–Verlet method, is an explicit
method for solving a system of second order differential equations. An
initial velocity vector :math:`V_0` is needed as well as the initial
condition :math:`X_0`.

.. math::
    f(t,X) = \ddot{X}

The first step is calculated with:

.. math::
    X_1 = X_0 + V_0 h + \frac{1}{2} f(X_0) h^2.

Subsequent steps are calculated with:

.. math::
    X_{n+1} = 2 X_{n} - X_{n-1} + f(t, X_n) h^2.

If subsequent velocities are needed, they can be calculated with:

.. math::
    V_n = \frac{X_{n+1} - X_{n-1}}{2 h}.
