from pathlib import Path

with open(Path(__file__).with_name("input_1_2.txt")) as file:
    rows = [[float(number) for number in line.split()] for line in file]

a = []
b = []
c = []
d = []

for row in rows:
    a.append(row[0])
    b.append(row[1])
    c.append(row[2])
    d.append(row[3])

n = len(b)

alpha = [0.0] * n
beta = [0.0] * n

alpha[0] = -c[0] / b[0]
beta[0] = d[0] / b[0]

for i in range(1, n):
    denominator = b[i] + a[i] * alpha[i - 1]
    alpha[i] = -c[i] / denominator
    beta[i] = (d[i] - a[i] * beta[i - 1]) / denominator

x = [0.0] * n
x[n - 1] = beta[n - 1]

for i in range(n - 2, -1, -1):
    x[i] = alpha[i] * x[i + 1] + beta[i]

print("Прогоночные коэффициенты:")
for i in range(n):
    print(f"alpha[{i + 1}] = {alpha[i]:.4f}, beta[{i + 1}] = {beta[i]:.4f}")

print("\nРешение:")
for i in range(n):
    print(f"x[{i + 1}] = {x[i]:.4f}")
