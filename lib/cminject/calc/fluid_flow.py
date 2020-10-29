#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
# later version. 
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details. 
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Calculation code for models of forces exerted on particles by flowing fluids.
"""

import numba
import numpy as np
from numpy import pi
from scipy.constants import Boltzmann

from cminject.calc import erf_vec_1d


@numba.jit(nopython=True, error_model='numpy')
def slip_correction_4k(p, r_p, T, ss):
    # The Sutherland constant for helium is 79.4 at reference temperature of 273.0 K.
    # I took a reference pressure of 1 Pascal, in which the mean free path of helium is 0.01754.
    # see J. Aerosol Sci. 1976 Vol. 7. pp 381-387 by Klaus Willeke
    knudsen = 0.01754 * 1 / (p * r_p) * (T / 273.0) * ((1 + 79.4 / 273.0) / (1 + 79.4 / T))
    s = 1 + knudsen * (1.246 + (0.42 * np.exp(-0.87 / knudsen)))
    return s * ss


@numba.jit(nopython=True, error_model='numpy')
def slip_correction_hutchins(mu, p, r_p, T_f, m_f, ss):
    """
    Slip correction for a spherical particle in a fluid according to:
    D. K. Hutchins, M. H. Harper, R. L. Felder, Slip correction measurements for solid spherical particles by modulated
    dynamic light scattering, Aerosol Sci. Techn. 22 (2) (1995) 202–218. doi:10.1080/02786829408959741.

    :param mu: The dynamic viscosity of the fluid.
    :param p: The pressure of the fluid.
    :param r_p: The radius of the particle.
    :param T_f: The temperature of the fluid.
    :param m_f: The mass of a single particle of the fluid.
    :param ss: A multiplicative correction factor.
    :return: The slip correction factor according to Hutchins.
    """
    knudsen = mu / (p * r_p) * np.sqrt(pi * Boltzmann * T_f / (2 * m_f))
    s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))
    return s * ss


@numba.jit(nopython=True, error_model='numpy')
def a_stokes(delta_v, mu, r_p, m_p, Cc):
    """
    Acceleration for a Stokes flow at a single point, 6*pi*mu*r_p*delta_v / (m_p * Cc).

    :param delta_v: The velocity difference (v_fluid - v_particle) at that point.
    :param mu: The dynamic viscosity of the fluid.
    :param r_p: The radius of the particle.
    :param m_p: The mass of the particle.
    :param Cc: A slip correction factor.
    :return: The acceleration according to Stokes.
    """
    return 6*pi * mu * r_p * delta_v / (m_p * Cc)


@numba.jit(nopython=True, error_model='numpy')
def a_brown_stokes(mu, T_f, r_p, rho_p, Cc, dt, n):
    """
    Acceleration for Brownian motion in a Stokes flow, according to Li and Ahmadi.

    :param mu: The dynamic viscosity of the fluid.
    :param T_f: The temperature of the fluid.
    :param r_p: The radius of the particle.
    :param rho_p: The density of the particle.
    :param Cc: A slip correction factor.
    :param dt: The time step.
    :param n: The number of dimensions to generate an acceleration for.
    :return: The acceleration vector.
    """
    s0 = 216 * mu * Boltzmann * T_f / (pi**2 * (2 * r_p)**5 * rho_p**2 * Cc)
    return np.random.normal(0.0, 1.0, n) * np.sqrt(pi * s0 / dt)


@numba.jit(nopython=True, error_model='numpy')
def a_roth(p, delta_v, m_f, T_f, T_p, r_p, m_p):
    """
    Acceleration for a microscopic force for aerosol transport described in Roth2020_.

    .. _Roth2020: https://arxiv.org/abs/2006.10652

    :param p: The pressure of the fluid.
    :param delta_v: The velocity difference (v_fluid - v_particle).
    :param m_f: The mass of a single particle of the fluid.
    :param T_f: The temperature of the fluid.
    :param T_p: The temperature of the particle.
    :param r_p: The radius of the particle.
    :param m_p: The mass of the particle.
    :return: The acceleration vector.
    """
    h_f = m_f / (2 * Boltzmann * T_f)
    h_p = m_f / (2 * Boltzmann * T_p)

    f_spec = np.where(
        np.abs(delta_v) < 10,  # For <10, we use Epstein's formula, for >10 the one by Nils Roth
        16/3 * p * np.sqrt(pi * h_f) * r_p**2 * delta_v,  # Epstein model
        -p * np.sqrt(pi) * r_p**2 * (
                -2 * np.exp(-h_f * delta_v**2) * np.sqrt(h_f) * delta_v * (1 + 2*h_f*delta_v**2) +
                np.sqrt(pi) * (1 - 4*h_f*delta_v**2 - 4*h_f**2*delta_v**4) *
                erf_vec_1d(np.sqrt(h_f) * delta_v)
        ) / (2 * h_f * delta_v**2)
    )
    return (f_spec + 1.8/3 * p * pi**(3/2) * h_f/np.sqrt(h_p) * r_p**2 * delta_v) / m_p


@numba.jit(nopython=True, error_model='numpy')
def a_brown_roth(p, m_f, T_f, T_p, r_p, m_p, dt, n):
    """
    Returns a random (Brownian) acceleration with a spectral intensity as described in Roth2020_,
    depending on particle and fluid properties at a certain position.

    :param p: The pressure at the point.
    :param m_f: The mass of a single particle of the fluid.
    :param T_f: The temperature of the fluid.
    :param T_p: The temperature of the particle.
    :param r_p: The radius of the particle.
    :param m_p: The mass of the particle.
    :param dt: The numerical time-step.
    :param n: The number of spatial dimensions.
    :return: An (n,)-shaped np.array containing the acceleration components in each spatial dimension.
    """
    h = m_f / (2 * Boltzmann * T_f)
    s0 = (16/3 + 2/3*pi*np.sqrt(T_p / T_f)) * np.sqrt(pi / h) * p * m_f * r_p**2
    return np.random.normal(0.0, 1.0, n) * np.sqrt(s0 / dt) / m_p
