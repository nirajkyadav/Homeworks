import numpy as np

np.set_printoptions(precision=2, suppress=True, floatmode="fixed")

SCALE = 800_000.0
A = 5.0 / 12.0          # Integration limits [-A, +A]


def B(xi):
    return np.array([
        [0.24 * xi],
        [0.6 * xi - 0.2],
        [-0.24 * xi],
        [0.6 * xi + 0.2]
        ])


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


def stiffness(n):
    t_ref, w = gauss_lobatto(n)
    K = np.zeros((4, 4))
    J = A                       # Jacobian for mapping [-1, 1] to [-A, A]

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

for n in (3, 5, 7):
    K = stiffness(n)
    err = np.max(np.abs(K - K_analytical))
    print(f"\n Gauss-Lobatto, {n} points")
    print(K)
    print(f"\n max|K - K_analytical| = {err:.3e}")
    print (f'================================================')

print("\n Analytical (exact integration)")
print(K_analytical)