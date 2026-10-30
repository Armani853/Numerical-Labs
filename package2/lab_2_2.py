import math


def read_input(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    n = int(lines[0])
    eps = float(lines[1])
    q = float(lines[2])
    x0 = list(map(float, lines[3].split()))
    return n, eps, q, x0



def f_vec(x):
    f1 = (x[0] ** 2) / 9.0 + (x[1] ** 2) / 2.25 - 1.0
    f2 = 3.0 * x[1] - math.exp(x[0]) - x[0]
    return [f1, f2]


def jacobian(x):
    j11 = (2.0 * x[0]) / 9.0
    j12 = (8.0 * x[1]) / 9.0
    j21 = -math.exp(x[0]) - 1.0
    j22 = 3.0
    return [[j11, j12], [j21, j22]]


def phi_vec(x):
    val = 2.25 * (1.0 - (x[0] ** 2) / 9.0)
    if val < 0:
        val = 0
    phi2 = math.sqrt(val)

    arg = 3.0 * x[1] - x[0]
    if arg <= 0:
        arg = 1e-5
    phi1 = math.log(arg)

    return [phi1, phi2]


def solve_gauss(A, b):
    n = len(b)
    M = []
    for i in range(n):
        M.append(A[i][:] + [b[i]])

    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(M[k][i]) > abs(M[max_row][i]):
                max_row = k
        M[i], M[max_row] = M[max_row], M[i]

        pivot = M[i][i]
        for j in range(i, n + 1):
            M[i][j] /= pivot

        for k in range(i + 1, n):
            factor = M[k][i]
            for j in range(i, n + 1):
                M[k][j] -= factor * M[i][j]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = M[i][n]
        for j in range(i + 1, n):
            x[i] -= M[i][j] * x[j]
    return x


def newton_method(x0, eps):
    print("НЬЮТОН:")
    x_curr = x0[:]
    k = 0
    print(f"k={k}: x1={x_curr[0]:.6f}, x2={x_curr[1]:.6f}")

    while True:
        k += 1
        F = f_vec(x_curr)
        J = jacobian(x_curr)

        neg_F = [-val for val in F]
        delta = solve_gauss(J, neg_F)

        x_next = [x_curr[0] + delta[0], x_curr[1] + delta[1]]
        diff = max(abs(x_next[0] - x_curr[0]), abs(x_next[1] - x_curr[1]))

        print(f"k={k}: x1={x_next[0]:.6f}, x2={x_next[1]:.6f}, diff={diff:.6e}")

        if diff < eps:
            return x_next, k

        x_curr = x_next


def simple_iteration(x0, eps, q):
    print("\nМПИ:")
    stop = ((1.0 - q) / q) * eps
    x_curr = x0[:]
    k = 0
    print(f"k={k}: x1={x_curr[0]:.6f}, x2={x_curr[1]:.6f}")

    while True:
        k += 1
        x_next = phi_vec(x_curr)
        diff = max(abs(x_next[0] - x_curr[0]), abs(x_next[1] - x_curr[1]))

        print(f"k={k}: x1={x_next[0]:.6f}, x2={x_next[1]:.6f}, diff={diff:.6e}")

        if diff <= stop:
            return x_next, k

        x_curr = x_next


def main():
    n, eps, q, x0 = read_input("input_2_2.txt")

    sol_n, k_n = newton_method(x0, eps)
    sol_i, k_i = simple_iteration(x0, eps, q)

    print("\nИТОГ:")
    print(f"Ньютон: x1={sol_n[0]:.6f}, x2={sol_n[1]:.6f}, итераций={k_n}")
    print(f"МПИ: x1={sol_i[0]:.6f}, x2={sol_i[1]:.6f}, итераций={k_i}")


if __name__ == "__main__":
    main()

