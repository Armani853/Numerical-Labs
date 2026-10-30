import math

def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    l = len(B)
    res = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(l):
                res[i][j] += A[i][k] * B[k][j]
    return res


def eye(n):
    E = [[0.0] * n for _ in range(n)]
    for i in range(n):
        E[i][i] = 1.0
    return E


def sign(x):
    if x >= 0:
        return 1.0
    return -1.0


def qr(A):
    n = len(A)
    R = [row[:] for row in A]
    Q = eye(n)

    for k in range(n - 1):
        b = [R[i][k] for i in range(k, n)]
        norm_b = math.sqrt(sum(x ** 2 for x in b))
        if norm_b == 0:
            continue

        v = [0.0] * len(b)
        v[0] = b[0] + sign(b[0]) * norm_b
        for i in range(1, len(b)):
            v[i] = b[i]

        v_sum = sum(x ** 2 for x in v)
        if v_sum == 0:
            continue

        m = len(b)
        H_sub = eye(m)
        for i in range(m):
            for j in range(m):
                H_sub[i][j] -= (2.0 * v[i] * v[j]) / v_sum

        H = eye(n)
        for i in range(m):
            for j in range(m):
                H[k + i][k + j] = H_sub[i][j]

        R = mat_mul(H, R)
        Q = mat_mul(Q, H)

    return Q, R


def solve_2x2(a11, a12, a21, a22):
    tr = a11 + a22
    det = a11 * a22 - a12 * a21
    d = tr ** 2 - 4 * det
    if d >= 0:
        return (tr + math.sqrt(d)) / 2.0, (tr - math.sqrt(d)) / 2.0
    else:
        real = tr / 2.0
        imag = math.sqrt(-d) / 2.0
        return complex(real, imag), complex(real, -imag)


def get_eigenvalues(A, eps=0.01):
    n = len(A)
    A_k = [row[:] for row in A]

    for _ in range(500):
        Q, R = qr(A_k)
        A_next = mat_mul(R, Q)

        diff = 0.0
        for r in range(n):
            for c in range(n):
                diff += (A_next[r][c] - A_k[r][c]) ** 2
        diff = math.sqrt(diff)

        A_k = A_next

        if diff < eps:
            break

    res = []
    i = 0
    while i < n:
        if i == n - 1:
            res.append(A_k[i][i])
            i += 1
        elif abs(A_k[i + 1][i]) < eps:
            res.append(A_k[i][i])
            i += 1
        else:
            l1, l2 = solve_2x2(A_k[i][i], A_k[i][i + 1], A_k[i + 1][i], A_k[i + 1][i + 1])
            res.append(l1)
            res.append(l2)
            i += 2

    return res, A_k


def print_matrix(M):
    for row in M:
        print(" ".join(f"{x:10.4f}" for x in row))


f = open("input_1_5.txt", "r")
lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
f.close()

eps = float(lines[0])
A = []
for line in lines[1:]:
    A.append([float(x) for x in line.split()])

print("Eps:", eps)
print("Matrix A:")
print_matrix(A)

Q, R = qr(A)
print("\nMatrix Q:")
print_matrix(Q)

print("\nMatrix R:")
print_matrix(R)

QR = mat_mul(Q, R)
print("\nCheck Q * R:")
print_matrix(QR)

vals, A_final = get_eigenvalues(A, eps)

print("\nEigenvalues:")
for i, val in enumerate(vals):
    if isinstance(val, complex):
        print(f"lambda_{i+1} = {val.real:.4f} + {val.imag:.4f}j")
    else:
        print(f"lambda_{i+1} = {val:.4f}")

