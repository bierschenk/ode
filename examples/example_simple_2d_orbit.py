
import ode
import math
import matplotlib.pyplot as plt


def dx_orbit_sys(t, X):
    m_earth = 5.97237*(10**24)  # kg
    m_moon = 7.342*(10**22)  # kg
    G = 6.67408*(10**-11)  # m**3 kg**−1 s**−2
    distance = ((X[0] - X[2])**2 + (X[1] - X[3])**2)**(1/2)
    force = G*m_earth*m_moon/(distance**2)
    direction = math.atan2(X[3] - X[1], X[2] - X[0])
    a_earth = force/m_earth
    a_moon = force/m_moon
    dX = [
            X[4],  # Vx Earth
            X[5],  # Vy Earth
            X[6],  # Vx Moon
            X[7],  # Vy Moon
            a_earth*math.cos(direction),
            a_earth*math.sin(direction),
            a_moon*math.cos(direction+math.pi),
            a_moon*math.sin(direction+math.pi),
            ]
    return dX


def ddx_orbit_sys(t, X):
    m_earth = 5.97237*(10**24)  # kg
    m_moon = 7.342*(10**22)  # kg
    G = 6.67408*(10**-11)  # m**3 kg**−1 s**−2
    distance = ((X[0] - X[2])**2 + (X[1] - X[3])**2)**(1/2)
    force = G*m_earth*m_moon/(distance**2)
    direction = math.atan2(X[3] - X[1], X[2] - X[0])
    a_earth = force/m_earth
    a_moon = force/m_moon
    ddX = [
            a_earth*math.cos(direction),
            a_earth*math.sin(direction),
            a_moon*math.cos(direction+math.pi),
            a_moon*math.sin(direction+math.pi),
            ]
    return ddX


t_range = [0, 10000000]
t_step = 1000

X_0 = [
        0,          # Px⋅Earth (m)
        0,          # Py Earth (m)
        384000000,  # Px Moon (m)
        0,          # Py Moon (m)
        0,          # Vx Earth (m/s)
        -12.5637,   # Vy Earth (m/s)
        0,          # Vx Moon (m/s)
        1022,       # Vy Moon (m/s)
        ]


t_euler, x_euler = ode.euler(
        dfun=dx_orbit_sys,
        xzero=X_0,
        timerange=t_range,
        timestep=t_step,
        )

t_backward_euler, x_backward_euler = ode.backwardeuler(
        dfun=dx_orbit_sys,
        xzero=X_0,
        timerange=t_range,
        timestep=t_step,
        )

e1,  e2,  e3,  e4,  e5,  e6,  e7,  e8 = x_euler
be1, be2, be3, be4, be5, be6, be7, be8 = x_backward_euler

x_verlet = [
        0,          # Px⋅Earth (m)
        0,          # Py Earth (m)
        384000000,  # Px Moon (m)
        0,          # Py Moon (m)
        ]

v_verlet_0 = [
        0,          # Vx Earth (m/s)
        -12.5637,   # Vy Earth (m/s)
        0,          # Vx Moon (m/s)
        1022,       # Vy Moon (m/s)
        ]

t_verlet, x_verlet, v_verlet = ode.verlet(
        dfun=ddx_orbit_sys, xzero=x_verlet, vzero=v_verlet_0,
        timerange=t_range, timestep=t_step)

v1, v2, v3, v4 = x_verlet

plt.plot(e1, e2, label='Earth, Euler')
plt.plot(e3, e4, label='Moon, Euler')
plt.plot(be1, be2, label='Earth, Backward Euler')
plt.plot(be3, be4, label='Moon, Backward Euler')
plt.plot(v1, v2, label='Earth, Verlet')
plt.plot(v3, v4, label='Moon, Verlet')
plt.legend()
plt.show()
