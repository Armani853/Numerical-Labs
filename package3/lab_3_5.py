import math
from pathlib import Path

def f(x):
    return 1.0 / (x**4 + 16.0)

def rect_left(a, b, h):
    n = int((b - a) / h + 0.5)
    s = 0.0
    for i in range(n):
        s += f(a + i * h)
    return s * h

def rect_right(a, b, h):
    n = int((b - a) / h + 0.5)
    s = 0.0
    for i in range(1, n + 1):
        s += f(a + i * h)
    return s * h

def rect_mid(a, b, h):
    n = int((b - a) / h + 0.5)
    s = 0.0
    for i in range(n):
        s += f(a + i * h + h / 2.0)
    return s * h

def trapezoid(a, b, h):
    n = int((b - a) / h + 0.5)
    s = (f(a) + f(b)) / 2.0
    for i in range(1, n):
        s += f(a + i * h)
    return s * h

def simpson(a, b, h):
    n = int((b - a) / h + 0.5)
    s = f(a) + f(b)
    for i in range(1, n):
        if i % 2 == 1:
            s += 4.0 * f(a + i * h)
        else:
            s += 2.0 * f(a + i * h)
    return s * h / 3.0

f_in = open(Path(__file__).with_name("input_3_5.txt"), "r")
data = f_in.read().split()
f_in.close()

a = float(data[0])
b = float(data[1])
h1 = float(data[2])
h2 = float(data[3])

methods = [
    ("Левые прямоугольники", rect_left, 1),
    ("Правые прямоугольники", rect_right, 1),
    ("Средние прямоугольники", rect_mid, 2),
    ("Метод трапеций", trapezoid, 2),
    ("Метод Симпсона", simpson, 4)
]

print("Результаты вычислений:")
print(f"{'Метод':<23} | {'I(h1)':<8} | {'I(h2)':<8} | {'p':<1} | {'I(Рунге)':<8} | {'Погрешность':<8}")
print("-" * 71)

for name, algo, p in methods:
    i1 = algo(a, b, h1)
    i2 = algo(a, b, h2)
    err = (i2 - i1) / (2.0**p - 1.0)
    i_runge = i2 + err
    print(f"{name:<23} | {i1:.5f} | {i2:.5f} | {p:<1} | {i_runge:.5f} | {abs(err):.5f}")

