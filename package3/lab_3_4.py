f_in = open("input_3_4.txt", "r")
lines = f_in.readlines()
f_in.close()

x_table = list(map(float, lines[0].split()))
y_table = list(map(float, lines[1].split()))
x_target = float(lines[2].strip())

idx = -1
for i in range(len(x_table)):
    if abs(x_table[i] - x_target) < 1e-9:
        idx = i
        break

if idx == -1:
    print("Точка не найдена в таблице")
elif idx == 0 or idx == len(x_table) - 1:
    print("Невозможно вычислить производные в крайней точке")
else:
    h = x_table[idx + 1] - x_table[idx]

    used_x = [x_table[idx - 1], x_table[idx], x_table[idx + 1]]
    used_y = [y_table[idx - 1], y_table[idx], y_table[idx + 1]]

    y_minus = y_table[idx - 1]
    y_curr = y_table[idx]
    y_plus = y_table[idx + 1]

    d1_left = (y_curr - y_minus) / h
    d1_right = (y_plus - y_curr) / h
    d1_central = (y_plus - y_minus) / (2.0 * h)

    d2 = (y_plus - 2.0 * y_curr + y_minus) / (h ** 2)

    print("Исходная таблица значений функции:")
    print("X:", x_table)
    print("Y:", y_table)
    print()
    print("Точка, в которой выполняется вычисление:", x_target)
    print()
    print("Узлы таблицы, использованные при вычислении:")
    print("X_nodes:", used_x)
    print("Y_nodes:", used_y)
    print()
    print("Значение первой производной:")
    print("  Левосторонняя:", round(d1_left, 4))
    print("  Правосторонняя:", round(d1_right, 4))
    print("  Центральная:", round(d1_central, 4))
    print()
    print("Значение второй производной:", round(d2, 4))

