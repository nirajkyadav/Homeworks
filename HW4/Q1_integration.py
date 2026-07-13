import numpy as np
from scipy.integrate import quad

J = 2.5            # Jacobian = L/2
L = 5.0            # Element length (L = 2 * J)

xi_min = -5 / 12
xi_max = 5 / 12

def get_kappa(xi):
    return -0.0024 * xi

def get_M(kappa):
    abs_k = abs(kappa)
    inside_sqrt = 2000 * abs_k - 1000_000 * kappa**2

    return np.sign(kappa) * 320 * np.sqrt(max(0.0, inside_sqrt))    # Avoid negative values inside square root

def get_dM_dk(kappa):
    abs_k = abs(kappa)
    num = 320 * (1000 - 1000_000 * abs_k)
    den = np.sqrt(2000 * abs_k - 1000_000 * kappa**2)

    return num / den if den != 0 else 0.0

def B_matrix(xi):
    return np.array([
        0.24 * xi,
        0.6 * xi - 0.2,
        -0.24 * xi,
        0.6 * xi + 0.2])

# Stiffness and Force Matrix
K2 = np.zeros((4, 4))
P2 = np.zeros(4)

for i in range(4):
    # P2 force components
    P2_integrand = lambda xi: B_matrix(xi)[i] * get_M(get_kappa(xi)) * J
    P2[i], _ = quad(P2_integrand, xi_min, xi_max)
    
    # K2 stiffness components
    for j in range(4):
        K2_integrand = lambda xi: B_matrix(xi)[i] * B_matrix(xi)[j] * get_dM_dk(get_kappa(xi)) * J
        K2[i, j], _ = quad(K2_integrand, xi_min, xi_max)

print("Stiffness Matrix Contribution (K2):")
print(np.round(K2, 2))

print("\nInternal Force Vector Contribution (P2):")
print(np.round(P2, 2))
