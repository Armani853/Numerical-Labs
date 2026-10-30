import math

def eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def transpose(A):
    n = len(A)
    m = len(A[0])
    res = [[0.0] * n for _ in range(m)]
    for i in range(n):
        for j in range(m):
            res[j][i] = A[i][j]
    return res

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

def print_mat(A, name):
    print(f" {name}")
    for row in A:
        print(" ".join(f"{x:10.4f}" for x in row))
    print()

def print_vec(v, name):
    print(f"{name}")
    print(" ".join(f"{x:10.4f}" for x in v))
    print()

def off_diag_norm(A):
    n = len(A)
    s = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            s += A[i][j] ** 2
    return math.sqrt(s)

def max_off_diag(A):
    n = len(A)
    max_val = -1.0
    p, q = 0, 1
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j]) > max_val:
                max_val = abs(A[i][j])
                p, q = i, j
    return p, q

def jacobi(A, eps):
    n = len(A)
    A_k = [row[:] for row in A]
    V = eye(n)

    while off_diag_norm(A_k) > eps:
        p, q = max_off_diag(A_k)

        if abs(A_k[p][p] - A_k[q][q]) < 1e-12:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2.0 * A_k[p][q] / (A_k[p][p] - A_k[q][q]))

        c = math.cos(phi)
        s = math.sin(phi)

        U = eye(n)
        U[p][p] = c
        U[q][q] = c
        U[p][q] = -s
        U[q][p] = s

        U_T = transpose(U)
        A_k = mat_mult(U_T, mat_mult(A_k, U))
        V = mat_mult(V, U)

    eigenvalues = [A_k[i][i] for i in range(n)]
    return eigenvalues, V

def main():
    f = open("input_1_4.txt", "r")
    lines = [line.strip() for line in f if line.strip()]
    f.close()

    n = int(lines[0])
    eps = float(lines[1])

    A = []
    for i in range(2, n + 2):
        row = list(map(float, lines[i].split()))
        A.append(row)

    print(f"eps = {eps}\n")
    print_mat(A, "A")

    eigenvalues, V = jacobi(A, eps)

    print_vec(eigenvalues, "Eigenvalues")
    print_mat(V, "V")

    Lambda = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Lambda[i][i] = eigenvalues[i]

    AV = mat_mult(A, V)
    VL = mat_mult(V, Lambda)

    print_mat(AV, "A * V")
    print_mat(VL, "V * Lambda")

if __name__ == "__main__":
    main()

