"""
==============================================================================
CEE 721: Nonlinear Structural Analysis

Date: July 13, 2026
Author: Niraj Yadav

HW 4, Problem 1: Internal Force and Tangent Stiffness Matrix

Units: kips, inches, ksi

==============================================================================
"""

import numpy as np

np.set_printoptions(precision=2, suppress=True, floatmode="fixed")

def B(xi):
    return np.array([
        [0.24 * xi],
        [0.6 * xi - 0.2],
        [-0.24 * xi],
        [0.6 * xi + 0.2]
    ])

def gauss_legendre(n):
    if n == 3:
        x = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])     # Gauss-Quadrature position
        w = np.array([5 / 9, 8 / 9, 5 / 9])                      # Gauss-Quadrature weights
    elif n == 5:
        a = (1 / 3) * np.sqrt(5 - 2 * np.sqrt(10 / 7))
        b = (1 / 3) * np.sqrt(5 + 2 * np.sqrt(10 / 7))
        x = np.array([-b, -a, 0.0, a, b])
        w = np.array([(322 - 13 * np.sqrt(70)) / 900,
                      (322 + 13 * np.sqrt(70)) / 900,
                      128 / 225,
                      (322 + 13 * np.sqrt(70)) / 900,
                      (322 - 13 * np.sqrt(70)) / 900])
    else:
        print("Only n = 3, 5 are provided")
        return None, None
    return x, w

def gauss_lobatto(n):
    if n == 3:
        x = np.array([-1.0, 0.0, 1.0])
        w = np.array([1 / 3, 4 / 3, 1 / 3])
    elif n == 5:
        x = np.array([-1.0, -np.sqrt(3 / 7), 0.0, np.sqrt(3 / 7), 1.0])
        w = np.array([1 / 10, 49 / 90, 32 / 45, 49 / 90, 1 / 10])
    elif n == 7:
        p1 = np.sqrt(5 / 11 - (2 / 11) * np.sqrt(5 / 3))   
        p2 = np.sqrt(5 / 11 + (2 / 11) * np.sqrt(5 / 3))
        x = np.array([-1.0, -p2, -p1, 0.0, p1, p2, 1.0])
        w = np.array([1 / 21,
                      (124 - 7 * np.sqrt(15)) / 350,
                      (124 + 7 * np.sqrt(15)) / 350,
                      256 / 525,
                      (124 + 7 * np.sqrt(15)) / 350,
                      (124 - 7 * np.sqrt(15)) / 350,
                      1 / 21])
    else:
        print("Only n = 3, 5, 7 are provided")
        return None, None
    return x, w

def constitutive(xi):
    kappa = -0.0024 * xi
    abs_kappa = abs(kappa)
    
    if abs_kappa > 0.001:
        M = np.sign(kappa) * 320.0
        dM_dk = 0.0
    else:
        term = 2000.0 * abs_kappa - 1000_000.0 * (abs_kappa**2)
        if term <= 0:
            M = 0.0
            dM_dk = 320_000.0 # Limit slope as kappa approaches 0
        else:
            M = np.sign(kappa) * 320.0 * np.sqrt(term)
            dM_dk = 320.0 * (1000.0 - 1000_000.0 * abs_kappa) / np.sqrt(term)

    return M, dM_dk

def numerical_integration(n, method):

    if method == 'quadrature':
        position, weights = gauss_legendre(n)
    else:
        position, weights = gauss_lobatto(n)

    P = np.zeros((4, 1))
    K = np.zeros((4, 4))
    J = 2.5                    # dx/dxi = Jacobian
    
    for xi, w in zip(position, weights):
        M, dM_dk = constitutive(xi)
        b = B(xi)
        P += w * b * M * J
        K += w * (b @ b.T) * dM_dk * J
    return P, K

# Internal Force and Tangent Stiffness Matrix
for scheme in [(3, 'quadrature'), (5, 'quadrature'), (3, 'lobatto'), (5, 'lobatto'), (7, 'lobatto')]:
    print("-------------------------------------------")
    print(f"{scheme[1]} ({scheme[0]} Points)")

    P, K = numerical_integration(scheme[0], scheme[1])
    
    print("Nodal Force Vector P:\n", P)
    print("Tangent Stiffness Matrix K_T:\n", K)