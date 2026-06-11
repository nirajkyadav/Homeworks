import numpy as np

np.set_printoptions(precision=2, suppress=True, floatmode="fixed")

SCALE = 800_000.0
A = 5.0 / 12.0          # integration limits [-A, +A]


def B(xi):
    return np.array([
        [0.24 * xi],
        [0.6 * xi - 0.2],
        [-0.24 * xi],
        [0.6 * xi + 0.2]
        ])


def gauss_legendre(n):
    if n == 3:
        x = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])     # Gauss-Quadrature nodes
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


def stiffness(n):
    t_ref, w = gauss_legendre(n)
    K = np.zeros((4, 4))
    J = A                             # Jacobian for mapping [-1, 1] to [-A, A]

    for t, wi in zip(t_ref, w):
        xi = A * t
        b = B(xi)
        b_bt = b @ b.T
        K = K + wi * b_bt * J

    return SCALE * K


K_analytical = np.array([
    [  2222.22,   5555.56,  -2222.22,   5555.56],
    [  5555.56,  40555.56,  -5555.56, -12777.78],
    [ -2222.22,  -5555.56,   2222.22,  -5555.56],
    [  5555.56, -12777.78,  -5555.56,  40555.56],
])


for n in (3, 5):
    K = stiffness(n)
    err = np.max(np.abs(K - K_analytical))
    print(f"\n Gauss-Quadrature, {n} points")
    print(K)
    print(f"\n max|K - K_analytical| = {err:.3e}")
    print (f'================================================')

print("\n Analytical (exact integration)")
print(K_analytical)
