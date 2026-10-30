def read_file(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()

    clean_lines = []
    for line in lines:
        if line.strip():
            clean_lines.append(line.strip())

    eps = float(clean_lines[0])

    matrix = []
    for line in clean_lines[1:]:
        row = []
        for val in line.split():
            row.append(float(val))
        matrix.append(row)

    n = len(matrix)
    A = []
    b = []
    for i in range(n):
        A.append(matrix[i][:n])
        b.append(matrix[i][n])

    return eps, A, b, n


def get_alpha_beta(A, b, n):
    alpha = []
    for i in range(n):
        alpha.append([0.0] * n)

    beta = [0.0] * n

    for i in range(n):
        beta[i] = b[i] / A[i][i]
        for j in range(n):
            if i != j:
                alpha[i][j] = -A[i][j] / A[i][i]

    max_sum = 0.0
    for i in range(n):
        row_sum = 0.0
        for j in range(n):
            row_sum += abs(alpha[i][j])
        if row_sum > max_sum:
            max_sum = row_sum

    return alpha, beta, max_sum


def mpi(alpha, beta, norm_a, n, eps):
    x_old = list(beta)
    k = 0
    coeff = norm_a / (1.0 - norm_a)

    while True:
        x_new = [0.0] * n
        k += 1

        for i in range(n):
            s = 0.0
            for j in range(n):
                s += alpha[i][j] * x_old[j]
            x_new[i] = beta[i] + s

        max_diff = 0.0
        for i in range(n):
            d = abs(x_new[i] - x_old[i])
            if d > max_diff:
                max_diff = d

        err = coeff * max_diff
        if err <= eps:
            return x_new, k

        x_old = x_new


def seidel(alpha, beta, norm_a, n, eps):
    x = list(beta)
    k = 0
    coeff = norm_a / (1.0 - norm_a)

    while True:
        x_old = list(x)
        k += 1

        for i in range(n):
            s = 0.0
            for j in range(n):
                s += alpha[i][j] * x[j]
            x[i] = beta[i] + s

        max_diff = 0.0
        for i in range(n):
            d = abs(x[i] - x_old[i])
            if d > max_diff:
                max_diff = d

        err = coeff * max_diff
        if err <= eps:
            return x, k


eps, A, b, n = read_file("input_1_3.txt")
alpha, beta, norm_a = get_alpha_beta(A, b, n)

x_mpi, k_mpi = mpi(alpha, beta, norm_a, n, eps)
x_seidel, k_seidel = seidel(alpha, beta, norm_a, n, eps)

print("Norma alpha =", round(norm_a, 4))
print("Eps =", eps)
print("-" * 30)

print("МПИ:")
print("Итераций:", k_mpi)
print("x =", [round(val, 4) for val in x_mpi])

print("\nМетод Зейделя:")
print("Итераций:", k_seidel)
print("x =", [round(val, 4) for val in x_seidel])

