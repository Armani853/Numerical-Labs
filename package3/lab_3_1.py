import math
import matplotlib.pyplot as plt

def f(x):
    return math.tan(x) + x

def lagrange(x, xs, ys):
    n = len(xs)
    res = 0.0
    for i in range(n):
        num = 1.0
        den = 1.0
        for j in range(n):
            if i != j:
                num *= (x - xs[j])
                den *= (xs[i] - xs[j])
        res += ys[i] * (num / den)
    return res

def lagrange_formula_str(xs, ys):
    n = len(xs)
    parts = []
    for i in range(n):
        if abs(ys[i]) < 1e-12:
            continue
        den = 1.0
        terms = []
        for j in range(n):
            if i != j:
                den *= (xs[i] - xs[j])
                terms.append(f"(x - {round(xs[j], 4)})")
        w = ys[i] / den
        parts.append(f"{w:+.6f} * " + " * ".join(terms))
    return "L(x) = " + "\n       ".join(parts)

def divided_diff_table(xs, ys):
    n = len(xs)
    table = [[0.0] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = ys[i]
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (xs[i + j] - xs[i])
    return table

def newton(x, xs, table):
    n = len(xs)
    res = table[0][0]
    prod = 1.0
    for i in range(1, n):
        prod *= (x - xs[i - 1])
        res += table[0][i] * prod
    return res

def newton_formula_str(xs, table):
    n = len(xs)
    s = f"P(x) = {table[0][0]:.6f}"
    terms = []
    for i in range(1, n):
        terms.append(f"(x - {round(xs[i - 1], 4)})")
        coef = table[0][i]
        s += f" {coef:+.6f} * " + " * ".join(terms)
    return s

def linspace(start, stop, num):
    step = (stop - start) / (num - 1)
    return [start + i * step for i in range(num)]


f_in = open("input_3_1.txt", "r")
lines = f_in.readlines()
f_in.close()

x_star = float(lines[0].strip())

all_variants = [
    [float(val) for val in lines[1].strip().split()],
    [float(val) for val in lines[2].strip().split()]
]
labels = ["А) Равномерная сетка", "Б) Неравномерная сетка"]

for step_idx, xs in enumerate(all_variants):
    ys = [f(x) for x in xs]
    n = len(xs)

    print(f"РАСЧЕТ ДЛЯ ПУНКТА {labels[step_idx]}")

    print("1. Исходные узлы интерполяции и значения функции:")
    print(" i        x_i             y_i")
    for i in range(n):
        print(f"{i:2d}  {xs[i]:14.8f}  {ys[i]:14.8f}")
    print()

    print("2. Интерполяционный многочлен Лагранжа:")
    print(lagrange_formula_str(xs, ys))
    print()

    table = divided_diff_table(xs, ys)
    print("3. Таблица разделённых разностей:")
    header = " i      x_i           y_i      " + "  ".join([f"Порядок {k}" for k in range(1, n)])
    print(header)
    for i in range(n):
        row_str = f"{i:2d}  {xs[i]:12.6f}"
        for j in range(n - i):
            row_str += f"  {table[i][j]:14.6f}"
        print(row_str)
    print()

    print("4. Интерполяционный многочлен Ньютона:")
    print(newton_formula_str(xs, table))
    print()

    y_lagr = lagrange(x_star, xs, ys)
    y_newt = newton(x_star, xs, table)
    y_exact = f(x_star)
    err_interp = abs(y_exact - y_newt)
    diff_methods = abs(y_lagr - y_newt)

    print(f"5. Значение многочлена Лагранжа в точке x*:  {y_lagr:.8f}")
    print(f"6. Значение многочлена Ньютона в точке x*:   {y_newt:.8f}")
    print(f"7. Точное значение функции f(x*):             {y_exact:.8f}")
    print(f"8. Абсолютная погрешность интерполяции:      {err_interp:.8e}")
    print(f"9. Разность между Лагранжем и Ньютоном:       {diff_methods:.8e}")
    print("\n")

    grid_x = linspace(xs[0], xs[-1], 200)
    grid_f = [f(x) for x in grid_x]
    grid_lagr = [lagrange(x, xs, ys) for x in grid_x]
    grid_newt = [newton(x, xs, table) for x in grid_x]

    plt.figure(figsize=(9, 5.5))
    plt.plot(grid_x, grid_f, label="f(x) = tg(x) + x", color="black", linewidth=2)
    plt.plot(grid_x, grid_lagr, label="Лагранж L(x)", color="blue", linestyle="--")
    plt.plot(grid_x, grid_newt, label="Ньютон P(x)", color="red", linestyle=":")
    plt.scatter(xs, ys, color="green", zorder=5, label="Узлы x_i")
    plt.scatter([x_star], [y_exact], color="black", marker="o", s=60, label="Точное f(x*)")
    plt.scatter([x_star], [y_newt], color="red", marker="x", s=80, label="Приближение (x*, P(x*))")
    plt.title(f"Интерполяция — Пункт {labels[step_idx]}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()

