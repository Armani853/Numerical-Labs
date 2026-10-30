import matplotlib.pyplot as plt
from pathlib import Path


def gauss(A, B):
    n = len(B)
    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > abs(A[max_row][i]):
                max_row = k
        A[i], A[max_row] = A[max_row], A[i]
        B[i], B[max_row] = B[max_row], B[i]

        for k in range(i + 1, n):
            factor = A[k][i] / A[i][i]
            for j in range(i, n):
                A[k][j] -= factor * A[i][j]
            B[k] -= factor * B[i]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = 0.0
        for j in range(i + 1, n):
            s += A[i][j] * x[j]
        x[i] = (B[i] - s) / A[i][i]
    return x


def make_system(x, y, m):
    n_eq = m + 1
    A = [[0.0] * n_eq for _ in range(n_eq)]
    B = [0.0] * n_eq
    for r in range(n_eq):
        for c in range(n_eq):
            power = r + c
            A[r][c] = sum(t ** power for t in x)
        B[r] = sum(y[i] * (x[i] ** r) for i in range(len(x)))
    return A, B


def eval_poly(coef, x_val):
    res = 0.0
    for i in range(len(coef)):
        res += coef[i] * (x_val ** i)
    return res


def calc_mse(coef, x, y):
    s = 0.0
    for i in range(len(x)):
        s += (eval_poly(coef, x[i]) - y[i]) ** 2
    return s / len(x)


def poly_str(coef):
    s = f"{coef[0]:.4f}"
    for i in range(1, len(coef)):
        s += f" {coef[i]:+.4f}*x^{i}"
    return s


f_in = open(Path(__file__).with_name("input_3_3.txt"), "r")
lines = f_in.readlines()
f_in.close()

x = list(map(float, lines[0].split()))
y = list(map(float, lines[1].split()))

A1, B1 = make_system(x, y, 1)
c1 = gauss(A1, B1)

A2, B2 = make_system(x, y, 2)
c2 = gauss(A2, B2)

A3, B3 = make_system(x, y, 3)
c3 = gauss(A3, B3)

mse1 = calc_mse(c1, x, y)
mse2 = calc_mse(c2, x, y)
mse3 = calc_mse(c3, x, y)

print("Коэффициенты многочлена 1-й степени:")
print(c1)
print("Коэффициенты многочлена 2-й степени:")
print(c2)
print("Коэффициенты多член 3-й степени:")
print(c3)
print()

print("Полученные многочлены в явном виде:")
print("P1(x) =", poly_str(c1))
print("P2(x) =", poly_str(c2))
print("P3(x) =", poly_str(c3))
print()

print(f"MSE1 = {mse1:.6f}")
print(f"MSE2 = {mse2:.6f}")
print(f"MSE3 = {mse3:.6f}")
print()

best = min(mse1, mse2, mse3)
if best == mse1:
    print("Лучше приближает многочлен 1-й степени")
elif best == mse2:
    print("Лучше приближает многочлен 2-й степени")
else:
    print("Лучше приближает многочлен 3-й степени")

x_min, x_max = min(x), max(x)
step = (x_max - x_min) / 199
xx = [x_min + i * step for i in range(200)]

yy1 = [eval_poly(c1, t) for t in xx]
yy2 = [eval_poly(c2, t) for t in xx]
yy3 = [eval_poly(c3, t) for t in xx]

plt.plot(xx, yy1, "--", label="1-я степень")
plt.plot(xx, yy2, "-.", label="2-я степень")
plt.plot(xx, yy3, "-", label="3-я степень")
plt.plot(x, y, "o", label="Точки")
plt.grid(True)
plt.legend()
plt.show()

