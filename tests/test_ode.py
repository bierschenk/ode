#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_ode
----------------------------------

Tests for `ode` module.
"""

import ode
from .test_data import *


# ----
# dot function definitions


def springmass(t, x):
    '''Given x_zero, return dx_zero for this system.
    First derivative. x is length 2.'''
    k = 5
    m = 1
    p, v = x
    dx_zero = [
           v,
           -p*k/m,
          ]
    return dx_zero


# -----
# TESTS


def test_euler():
    t_euler_raw, x_euler_raw = ode.euler(
        oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_euler_raw, v_euler_raw = zip(*x_euler_raw)
    t_euler = [round(x, 6) for x in t_euler_raw]
    p_euler = [round(x, 6) for x in p_euler_raw]
    v_euler = [round(x, 6) for x in v_euler_raw]
    t_ieuler_raw, x_ieuler_raw = zip(*list(ode.Euler(
        oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)))
    p_ieuler_raw, v_ieuler_raw = zip(*x_ieuler_raw)
    t_ieuler = [round(x, 6) for x in t_ieuler_raw]
    p_ieuler = [round(x, 6) for x in p_ieuler_raw]
    v_ieuler = [round(x, 6) for x in v_ieuler_raw]
    assert (t_euler, p_euler, v_euler) == (t_ieuler, p_ieuler, v_ieuler)


def test_Euler():
    t_test = [round(x, 6) for x in oscillator_euler_t]
    p_test = [round(x, 6) for x in oscillator_euler_x1]
    v_test = [round(x, 6) for x in oscillator_euler_x2]
    t_euler_raw, x_euler_raw = ode.euler(
            oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_euler_raw, v_euler_raw = zip(*x_euler_raw)
    t_euler = [round(x, 6) for x in t_euler_raw]
    p_euler = [round(x, 6) for x in p_euler_raw]
    v_euler = [round(x, 6) for x in v_euler_raw]
    assert (t_test, p_test, v_test) == (t_euler, p_euler, v_euler)


def test_ibackward_euler():
    t_back_e_raw, x_back_e_raw = ode.backwardeuler(
            oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_back_e_raw, v_back_e_raw = zip(*x_back_e_raw)
    t_back_e = [round(x, 6) for x in t_back_e_raw]
    p_back_e = [round(x, 6) for x in p_back_e_raw]
    v_back_e = [round(x, 6) for x in v_back_e_raw]
    t_iback_e_raw, x_iback_e_raw = zip(*list(ode.BackwardEuler(
            oscillator_1st_deriv, xzero=[0, 1],
            timerange=[0, 5], timestep=0.1)))
    p_iback_e_raw, v_iback_e_raw = zip(*x_iback_e_raw)
    t_iback_e = [round(x, 6) for x in t_iback_e_raw]
    p_iback_e = [round(x, 6) for x in p_iback_e_raw]
    v_iback_e = [round(x, 6) for x in v_iback_e_raw]
    assert (t_back_e, p_back_e, v_back_e) == (t_iback_e, p_iback_e, v_iback_e)


def test_backwardeuler():
    t_test_raw = [
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    t_test = [round(x, 8) for x in t_test_raw]
    p_test_raw = [
        1.0, 0.952380952382812, 0.861678004542233,
        0.734261958766686, 0.577948488565721, 0.401557160349537,
        0.214443649652305, 0.026028703767535, -0.154653563921079,
        -0.319367458676474, -0.461029860412492, -0.573992630619344,
        -0.654243238884226, -0.699517949668252, -0.709326343290485,
        -0.6848902256338, -0.629003912361704, -0.545826284849587,
        -0.440617768894858, -0.319437383754144, -0.188816189156835]
    p_test = [round(x, 8) for x in p_test_raw]
    v_test_raw = [
        0.0, -0.476190476171875, -0.907029478424391,
        -1.2741604577909, -1.56313470205942, -1.7639132822229,
        -1.87113510704121, -1.88414945892079, -1.80682267695975,
        -1.64713894762453, -1.41662401742452, -1.12962770212385,
        -0.80250608269295, -0.452747107871603, -0.098083936240023,
        0.244361176563023, 0.558863132730499, 0.831776275143007,
        1.05208515957978, 1.21180385144824, 1.30621194602042]
    v_test = [round(x, 8) for x in v_test_raw]
    t_back_e_raw, x_back_e_raw = ode.backwardeuler(
            dfun=springmass, xzero=[1, 0], timerange=[0, 2], timestep=0.1)
    p_back_e_raw, v_back_e_raw = zip(*x_back_e_raw)
    t_back_e = [round(x, 8) for x in t_back_e_raw]
    p_back_e = [round(x, 8) for x in p_back_e_raw]
    v_back_e = [round(x, 8) for x in v_back_e_raw]
    assert (t_test, p_test, v_test) == (t_back_e, p_back_e, v_back_e)


def test_Verlet():
    tv_raw, x_v_raw, v_v_raw = ode.verlet(
            dfun=oscillator_2nd_deriv, xzero=[0], vzero=[1],
            timerange=[0, 5], timestep=.1)
    xv_raw = list(zip(*x_v_raw))[0]
    vv_raw = list(zip(*v_v_raw))[0]
    tv = [round(x, 6) for x in tv_raw]
    xv = [round(x, 6) for x in xv_raw]
    vv = [round(x, 6) for x in vv_raw]
    tiv_raw, x_iv_raw, v_iv_raw = zip(*list(ode.Verlet(
            dfun=oscillator_2nd_deriv, xzero=[0], vzero=[1],
            timerange=[0, 5], timestep=.1)))
    xiv_raw = list(zip(*x_iv_raw))[0]
    viv_raw = list(zip(*v_iv_raw))[0]
    tiv = [round(x, 6) for x in tiv_raw]
    xiv = [round(x, 6) for x in xiv_raw]
    viv = [round(x, 6) for x in viv_raw]
    assert (tv, xv) == (tiv, xiv)


def test_verlet():
    t_test = [round(x, 6) for x in oscillator_verlet_t]
    x_test = [round(x, 6) for x in oscillator_verlet_x1]
    v_test = [round(x, 6) for x in oscillator_verlet_v1]
    t_v_raw, x_v_raw, v_v_raw = ode.verlet(
            dfun=oscillator_2nd_deriv, xzero=[0], vzero=[1],
            timerange=[0, 5], timestep=.1)
    xv_raw = list(zip(*x_v_raw))[0]
    vv_raw = list(zip(*v_v_raw))[0]
    t_v = [round(x, 6) for x in t_v_raw]
    x_v = [round(x, 6) for x in xv_raw]
    v_v = [round(x, 6) for x in vv_raw]
    assert t_test == t_v
    assert x_test == x_v
    assert v_test == v_v
