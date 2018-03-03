#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_ode
----------------------------------

Tests for `ode` module.
"""

import pytest
import ode



#-----
# dot function definitions

def spring_mass(*, X):
    'Given x_zero, return dx_zero for this system'
    k=5
    m=1
    p,v = X
    dx_zero = [
           v,
           -p*k/m,
          ]
    return dx_zero



#-----
# TESTS


def test_t_gen():
    t_test = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]
    t = ode._functions._t_gen(
            t_start = 0, t_end = 1.95, t_step = 0.1)
    t_rnd = [round(x,8) for x in t]
    assert t_rnd == t_test


def test_euler():
    t_test_raw = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
            1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    t_test = [round(x,4) for x in t_test_raw]
    p_test_raw = [1.0, 1.0, 0.95, 0.85, 0.7025, 0.5125, 0.287375,
            0.036625, -0.22849375, -0.49544375, -0.7509690625,
            -0.9817221875, -1.174926859375, -1.319045421875,
            -1.40441764140625, -1.42383758984375, -1.37303665621094,
            -1.25104384308594, -1.06039919715039, -0.807202359060547,
            -0.500985561113183]
    p_test = [round(x,8) for x in p_test_raw]
    v_test_raw = [0.0, -0.5, -1.0, -1.475, -1.9, -2.25125, -2.5075,
            -2.6511875, -2.6695, -2.555253125, -2.30753125,
            -1.93204671875, -1.441185625, -0.853722195312499,
            -0.194199484374999, 0.508009336328126, 1.21992813125,
            1.90644645935547, 2.53196838089844, 3.06216797947363,
            3.46576915900391]
    v_test = [round(x,8) for x in v_test_raw]
    t_euler_raw, x_euler_raw =  ode.euler(dot_func=spring_mass, x_zero=[1,0],
            t_range=[0,2], t_step=0.1)
    p_euler_raw, v_euler_raw = zip(*x_euler_raw)
    t_euler = [round(x,8) for x in t_euler_raw]
    p_euler = [round(x,8) for x in p_euler_raw]
    v_euler = [round(x,8) for x in v_euler_raw]
    assert (t_test, p_test, v_test) == (t_euler, p_euler, v_euler)


def test_backward_euler():
    t_test_raw = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    t_test = [round(x,4) for x in t_test_raw]
    p_test_raw = [1.0, 0.952380952382812, 0.861678004542233,
        0.734261958766686, 0.577948488565721, 0.401557160349537,
        0.214443649652305, 0.026028703767535, -0.154653563921079,
        -0.319367458676474, -0.461029860412492, -0.573992630619344,
        -0.654243238884226, -0.699517949668252, -0.709326343290485,
        -0.6848902256338, -0.629003912361704, -0.545826284849587,
        -0.440617768894858, -0.319437383754144, -0.188816189156835]
    p_test = [round(x,8) for x in p_test_raw]
    v_test_raw = [0.0, -0.476190476171875, -0.907029478424391,
        -1.2741604577909, -1.56313470205942, -1.7639132822229,
        -1.87113510704121, -1.88414945892079, -1.80682267695975,
        -1.64713894762453, -1.41662401742452, -1.12962770212385,
        -0.80250608269295, -0.452747107871603, -0.098083936240023,
        0.244361176563023, 0.558863132730499, 0.831776275143007,
        1.05208515957978, 1.21180385144824, 1.30621194602042]
    v_test = [round(x,8) for x in v_test_raw]
    t_back_e_raw, x_back_e_raw =  ode.backward_euler(dot_func=spring_mass, 
            x_zero=[1,0], t_range=[0,2], t_step=0.1)
    p_back_e_raw, v_back_e_raw = zip(*x_back_e_raw)
    t_back_e = [round(x,8) for x in t_back_e_raw]
    p_back_e = [round(x,8) for x in p_back_e_raw]
    v_back_e = [round(x,8) for x in v_back_e_raw]
    assert (t_test, p_test, v_test) == (t_back_e, p_back_e, v_back_e)




