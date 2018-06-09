
import ode
from math import atan2, cos, pi, sin, sqrt
import matplotlib.pyplot as plt


def dx_orbit_sys(t, X):
    '''X = [
    m1x, m1y,
    m2x, m2y,
    m3x, m3y,
    m4x, m4y,
    m1vx, m1vy,
    m2vx, m2vy,
    m3vx, m3vy,
    m4vx, m4vy
    ]
    '''
    (m1x, m1y,
     m2x, m2y,
     m3x, m3y,
     m4x, m4y,
     m1vx, m1vy,
     m2vx, m2vy,
     m3vx, m3vy,
     m4vx, m4vy) = X
    m_moon1 = 7.342*(10**22)  # kg
    m_moon2 = 7.342*(10**22)  # kg
    m_moon3 = 7.342*(10**22)  # kg
    m_moon4 = 7.342*(10**22)  # kg
    G = 6.67408*(10**-11)  # m**3 kg**−1 s**−2
    dm12 = sqrt((m1x - m2x)**2 + (m1y - m2y)**2)
    dm13 = sqrt((m1x - m3x)**2 + (m1y - m3y)**2)
    dm14 = sqrt((m1x - m4x)**2 + (m1y - m4y)**2)
    dm23 = sqrt((m2x - m3x)**2 + (m2y - m3y)**2)
    dm24 = sqrt((m2x - m4x)**2 + (m2y - m4y)**2)
    dm34 = sqrt((m3x - m4x)**2 + (m3y - m4y)**2)
    f12 = G * m_moon1 * m_moon2 / (dm12 * dm12)
    f13 = G * m_moon1 * m_moon3 / (dm13 * dm13)
    f14 = G * m_moon1 * m_moon4 / (dm14 * dm14)
    f23 = G * m_moon2 * m_moon3 / (dm23 * dm23)
    f24 = G * m_moon2 * m_moon4 / (dm24 * dm24)
    f34 = G * m_moon3 * m_moon4 / (dm34 * dm34)
    dr12 = atan2(m2y - m1y, m2x - m1x)
    dr13 = atan2(m3y - m1y, m3x - m1x)
    dr14 = atan2(m4y - m1y, m4x - m1x)
    dr23 = atan2(m3y - m2y, m3x - m2x)
    dr24 = atan2(m4y - m2y, m4x - m2x)
    dr34 = atan2(m4y - m3y, m4x - m3x)
    f1x = f12 * cos(dr12) + f13 * cos(dr13) + f14 * cos(dr14)
    f1y = f12 * sin(dr12) + f13 * sin(dr13) + f14 * sin(dr14)
    f2x = f12 * cos(dr12 + pi) + f23 * cos(dr23) + f24 * cos(dr24)
    f2y = f12 * sin(dr12 + pi) + f23 * sin(dr23) + f24 * sin(dr24)
    f3x = f13 * cos(dr13 + pi) + f23 * cos(dr23 + pi) + f34 * cos(dr34)
    f3y = f13 * sin(dr13 + pi) + f23 * sin(dr23 + pi) + f34 * sin(dr34)
    f4x = f14 * cos(dr14 + pi) + f24 * cos(dr24 + pi) + f34 * cos(dr34 + pi)
    f4y = f14 * sin(dr14 + pi) + f24 * sin(dr24 + pi) + f34 * sin(dr34 + pi)
    dX = [
            m1vx,
            m1vy,
            m2vx,
            m2vy,
            m3vx,
            m3vy,
            m4vx,
            m4vy,
            f1x / m_moon1,
            f1y / m_moon1,
            f2x / m_moon2,
            f2y / m_moon2,
            f3x / m_moon3,
            f3y / m_moon3,
            f4x / m_moon4,
            f4y / m_moon4,
            ]
    return dX


# def ddx_orbit_sys(*, X):
#     m_earth = 5.97237*(10**24)  # kg
#     m_moon = 7.342*(10**22)  # kg
#     G = 6.67408*(10**-11)  # m**3 kg**−1 s**−2
#     distance = ((X[0] - X[2])**2 + (X[1] - X[3])**2)**(1/2)
#     force = G*m_earth*m_moon/(distance**2)
#     direction = math.atan2(X[3] - X[1], X[2] - X[0])
#     a_earth = force/m_earth
#     a_moon = force/m_moon
#     ddX = [
#             a_earth*math.cos(direction),
#             a_earth*math.sin(direction),
#             a_moon*math.cos(direction+math.pi),
#             a_moon*math.sin(direction+math.pi),
#             ]
#     return ddX


t_range = [0, 100000000]
t_step = 10000

X_0 = [
        384000000,  # m1x
        0,          # m1y
        0,          # m2x
        384000000,  # m2y
        -384000000,  # m3x
        0,          # m3y
        0,          # m4x
        -384000000,  # m4y
        0,          # m1vx
        102.2,       # m1vy
        -102.2,      # m2vx
        0,          # m2vy
        0,          # m3vx
        -102.2,      # m3vy
        102.2,       # m4vx
        0,          # m4vy
        ]


t_euler, x_euler = ode.euler(
        dfun=dx_orbit_sys,
        xzero=X_0,
        timerange=t_range,
        timestep=t_step,
        )

# t_backward_euler, x_backward_euler = ode.backward_euler(
#         dot_func=dx_orbit_sys,
#         x_zero=X_0,
#         t_range=t_range,
#         t_step=t_step,
#         )

(em1x, em1y,
 em2x, em2y,
 em3x, em3y,
 em4x, em4y,
 em1vx, em1vy,
 em2vx, em2vy,
 em3vx, em3vy,
 em4vx, em4vy) = zip(*x_euler)

# (bem1x, bem1y,
#  bem2x, bem2y,
#  bem3x, bem3y,
#  bem4x, bem4y,
#  bem1vx, bem1vy,
#  bem2vx, bem2vy,
#  bem3vx, bem3vy,
#  bem4vx, bem4vy) = zip(*x_backward_euler)

# t_verlet, x_verlet = ode.verlet(
#         ddot_func=ddx_orbit_sys, x_zero=x_verlet, v_zero=v_verlet_0,
#         t_range=t_range, t_step=t_step)

# v1, v2, v3, v4 = zip(*x_verlet)

plt.plot(em1x, em1y, label='M1-E')
plt.plot(em2x, em2y, label='M2-E')
plt.plot(em3x, em3y, label='M3-E')
plt.plot(em4x, em4y, label='M4-E')
# plt.plot(bem1x, bem1y, label='M1-BE')
# plt.plot(bem2x, bem2y, label='M2-BE')
# plt.plot(bem3x, bem3y, label='M3-BE')
# plt.plot(bem4x, bem4y, label='M4-BE')
# plt.plot(v1, v2, label='Earth, Verlet')
# plt.plot(v3, v4, label='Moon, Verlet')
plt.legend()
plt.show()
