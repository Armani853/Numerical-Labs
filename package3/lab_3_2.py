
import matplotlib.pyplot as plt
from pathlib import Path


def progonka(a, b, c, d):
    n = len(d)
    p = [0.0] * n
    q = [0.0] * n
    p[0] = c[0] / b[0]
    q[0] = d[0] / b[0]

    for i in range(1, n):
        t = b[i] - a[i] * p[i - 1]
        if abs(t) < 1e-12:
            raise ValueError('Деление на ноль в прогонке')
        p[i] = c[i] / t
        q[i] = (d[i] - a[i] * q[i - 1]) / t

    x = [0.0] * n
    x[n - 1] = q[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = q[i] - p[i] * x[i + 1]
    return x


def spline_coeffs(x, y):
    n = len(x) - 1
    h = [x[i + 1] - x[i] for i in range(n)]
    m = n - 1
    lower = [0.0] * m
    diag = [0.0] * m
    upper = [0.0] * m
    rhs = [0.0] * m

    for i in range(m):
        k = i + 1
        diag[i] = 2 * (h[k - 1] + h[k])
        if i > 0:
            lower[i] = h[k - 1]
        if i < m - 1:
            upper[i] = h[k]
        rhs[i] = 3 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1])

    c = [0.0] + (progonka(lower, diag, upper, rhs) if m > 0 else []) + [0.0]
    a = []
    b = []
    d = []

    for i in range(n):
        a.append(y[i])
        b.append((y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3)
        d.append((c[i + 1] - c[i]) / (3 * h[i]))

    return a, b, c, d


def spline(t, x, a, b, c, d):
    n = len(x) - 1
    if t < x[0] or t > x[-1]:
        raise ValueError('Точка вне отрезка')

    for i in range(n):
        if x[i] <= t <= x[i + 1]:
            dx = t - x[i]
            s = a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
            return s, i + 1


f_in = open(Path(__file__).with_name('input_3_2.txt'), 'r')
lines = f_in.readlines()
f_in.close()

x_star = float(lines[0].strip())
x = list(map(float, lines[1].split()))
y = list(map(float, lines[2].split()))

if len(x) != len(y) or len(x) < 2:
    raise ValueError('Неверное количество узлов')
if any(x[i] >= x[i + 1] for i in range(len(x) - 1)):
    raise ValueError('Узлы должны идти по возрастанию')

a, b, c, d = spline_coeffs(x, y)
n = len(x) - 1

print('Исходные данные:')
print('i        x           y')
for i in range(n + 1):
    print(f'{i}  {x[i]:10.5f}  {y[i]:10.5f}')

print('\nКоэффициенты сплайна:')
print('i      интервал          a           b           c           d')
for i in range(n):
    print(f'{i+1}  [{x[i]:.2f}; {x[i+1]:.2f}]  {a[i]:11.7f} {b[i]:11.7f} {c[i]:11.7f} {d[i]:11.7f}')

s, interval = spline(x_star, x, a, b, c, d)
print('\nx* =', x_star)
print('Интервал:', interval)
print('S(x*) =', round(s, 10))

print('\nПроверка узлов:')
print('i       y          S(x)         ошибка')
for i in range(n + 1):
    value, _ = spline(x[i], x, a, b, c, d)
    print(f'{i}  {y[i]:11.7f} {value:11.7f} {abs(y[i] - value):12.3e}')

left = 2 * c[0]
right = 2 * c[-2] + 6 * d[-1] * (x[-1] - x[-2])
print('\nВторая производная слева:', left)
print('Вторая производная справа:', right)

xx = [x[0] + i * (x[-1] - x[0]) / 400 for i in range(401)]
yy = [spline(t, x, a, b, c, d)[0] for t in xx]

plt.plot(xx, yy, label='Кубический сплайн')
plt.plot(x, y, 'o', label='Узлы')
plt.scatter([x_star], [s], color='red', label='S(x*)')
plt.axvline(x_star, linestyle=':', color='gray')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.savefig(Path(__file__).with_name('graph_3_2.png'), dpi=200)
plt.show()

