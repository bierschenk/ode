#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_ode
----------------------------------

Tests for `ode` module.
"""

import ode
from .test_data import *


# -----
# TESTS


def test_euler():
    t_euler_raw, x_euler_raw = ode.euler(
        oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_euler_raw, v_euler_raw = x_euler_raw
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
    p_euler_raw, v_euler_raw = x_euler_raw
    t_euler = [round(x, 6) for x in t_euler_raw]
    p_euler = [round(x, 6) for x in p_euler_raw]
    v_euler = [round(x, 6) for x in v_euler_raw]
    assert (t_test, p_test, v_test) == (t_euler, p_euler, v_euler)


def test_backwardeuler():
    t_back_e_raw, x_back_e_raw = ode.backwardeuler(
            oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_back_e_raw, v_back_e_raw = x_back_e_raw
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


def test_BackwardEuler():
    t_test = [round(x, 6) for x in oscillator_backwardeuler_t]
    p_test = [round(x, 6) for x in oscillator_backwardeuler_x1]
    v_test = [round(x, 6) for x in oscillator_backwardeuler_x2]
    t_back_e_raw, x_back_e_raw = ode.backwardeuler(
            oscillator_1st_deriv, xzero=[0, 1], timerange=[0, 5], timestep=0.1)
    p_back_e_raw, v_back_e_raw = x_back_e_raw
    t_back_e = [round(x, 6) for x in t_back_e_raw]
    p_back_e = [round(x, 6) for x in p_back_e_raw]
    v_back_e = [round(x, 6) for x in v_back_e_raw]
    assert (t_test, p_test, v_test) == (t_back_e, p_back_e, v_back_e)


def test_Verlet():
    tv_raw, x_v_raw, v_v_raw = ode.verlet(
            dfun=oscillator_2nd_deriv, xzero=[0], vzero=[1],
            timerange=[0, 5], timestep=.1)
    xv_raw = list(x_v_raw)[0]
    vv_raw = list(v_v_raw)[0]
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
    xv_raw = list(x_v_raw)[0]
    vv_raw = list(v_v_raw)[0]
    t_v = [round(x, 6) for x in t_v_raw]
    x_v = [round(x, 6) for x in xv_raw]
    v_v = [round(x, 6) for x in vv_raw]
    assert t_test == t_v
    assert x_test == x_v
    assert v_test == v_v
