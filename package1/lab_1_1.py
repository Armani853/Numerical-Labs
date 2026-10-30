from pathlib import Path


def mat_mult(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for p in range(k):
                C[i][j] += A[i][p] * B[p][j]
    return C


def mat_vec_mult(A, x):
    n = len(A)
    res = [0.0] * n
    for i in range(n):
        for j in range(n):
            res[i] += A[i][j] * x[j]
    return res


def print_mat(A, name):
    print(f"{name}")
    for row in A:
        print(" ".join(f"{x:10.4f}" for x in row))
    print()


def print_vec(v, name):
    print(f"{name}")
    print(" ".join(f"{x:10.4f}" for x in v))
    print()


def lu_decomp(A):
    n = len(A)
    U = [row[:] for row in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    P = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    swaps = 0

    for k in range(n):
        max_i = k
        for i in range(k + 1, n):
            if abs(U[i][k]) > abs(U[max_i][k]):
                max_i = i

        if max_i != k:
            U[k], U[max_i] = U[max_i], U[k]
            P[k], P[max_i] = P[max_i], P[k]
            swaps += 1
            for j in range(k):
                L[k][j], L[max_i][j] = L[max_i][j], L[k][j]

        for i in range(k + 1, n):
            factor = U[i][k] / U[k][k]
            L[i][k] = factor
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] -= factor * U[k][j]

    return L, U, P, swaps


def solve_lu(L, U, P, b):
    n = len(L)
    Pb = mat_vec_mult(P, b)

    z = [0.0] * n
    for i in range(n):
        s = sum(L[i][j] * z[j] for j in range(i))
        z[i] = Pb[i] - s

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (z[i] - s) / U[i][i]

    return x


def get_det(U, swaps):
    n = len(U)
    det = 1.0 if swaps % 2 == 0 else -1.0
    for i in range(n):
        det *= U[i][i]
    return det


def get_inverse(L, U, P):
    n = len(L)
    A_inv = [[0.0] * n for _ in range(n)]
    for j in range(n):
        e = [1.0 if i == j else 0.0 for i in range(n)]
        col = solve_lu(L, U, P, e)
        for i in range(n):
            A_inv[i][j] = col[i]
    return A_inv


def main():
    f = open(Path(__file__).with_name("input_1_1.txt"), "r")
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines if line.strip()]
    n = int(lines[0])

    A = []
    b = []
    for i in range(1, n + 1):
        row = list(map(float, lines[i].split()))
        A.append(row[:n])
        b.append(row[n])

    print_mat(A, "A")
    print_vec(b, "b")

    L, U, P, swaps = lu_decomp(A)

    print_mat(L, "L")
    print_mat(U, "U")

    LU = mat_mult(L, U)
    print_mat(LU, "L * U")

    PA = mat_mult(P, A)
    print_mat(PA, "P * A")

    x = solve_lu(L, U, P, b)
    print_vec(x, "x")

    det = get_det(U, swaps)
    print(f"det(A) = {det:.4f}\n")

    A_inv = get_inverse(L, U, P)
    print_mat(A_inv, "A^-1")

    E_check = mat_mult(A, A_inv)
    print_mat(E_check, "A * A^-1")


if __name__ == "__main__":
    main()

