import math

def f(x):
    return x**3 - 2 * x**2 - 10 * x + 15

def df(x):
    return 3 * x**2 - 4 * x - 10

def d2f(x):
    return 6 * x - 4

def phi(x):
    return 0.1 * x**3 - 0.2 * x**2 + 1.5

def newton(a, b, eps):
    if f(a) * d2f(a) > 0:
        x0 = a
    elif f(b) * d2f(b) > 0:
        x0 = b
    else:
        x0 = (a + b) / 2.0

    res = []
    k = 0
    x = x0

    while True:
        y = f(x)
        dy = df(x)
        x_new = x - y / dy
        d = abs(x_new - x)

        res.append((k, x, y, d))

        if d < eps:
            res.append((k + 1, x_new, f(x_new), d))
            break

        x = x_new
        k += 1

    return x_new, res

def simple_iterations(a, b, eps, q):
    x = (a + b) / 2.0
    res = []
    k = 0
    st = ((1.0 - q) / q) * eps

    while True:
        x_new = phi(x)
        d = abs(x_new - x)
        y = f(x)

        res.append((k, x, y, d))

        if d < st:
            res.append((k + 1, x_new, f(x_new), d))
            break

        x = x_new
        k += 1

    return x_new, res

f_in = open("input_2_1.txt", "r")
lines = f_in.readlines()
f_in.close()

a = float(lines[0].strip())
b = float(lines[1].strip())
eps = float(lines[2].strip())
q = float(lines[3].strip())

r1, t1 = newton(a, b, eps)
print("Метод Ньютона")
print("k  x  f(x)  diff")
for r in t1:
    print(str(r[0]) + "  " + str(round(r[1], 6)) + "  " + str(r[2]) + "  " + str(r[3]))
print("Корень:", round(r1, 6))
print()

r2, t2 = simple_iterations(a, b, eps, q)
print("Метод простой итерации")
print("k  x  f(x)  diff")
for r in t2:
    print(str(r[0]) + "  " + str(round(r[1], 6)) + "  " + str(r[2]) + "  " + str(r[3]))
print("Корень:", round(r2, 6))

