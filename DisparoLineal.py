from tabulate import tabulate
import warnings
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)


def u_2_equation(x, u_1, u_2):
    # ! CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    return 4 - (1 / 4) * pow(x, 3) - u_1 * u_2


def runge_kutta_fourth_order(x, u_1, u_2, h):
    u_1_eval = u_2
    u_2_eval = u_2_equation(x, u_1, u_2)

    k_11 = h * u_1_eval
    k_12 = h * u_2_eval
    k_21 = h * (u_2 + k_12 / 2)
    # ! CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_22 = h * (4 - (1 / 4) * pow(x + h / 2, 3) - (u_1 + k_11 / 2) * (u_2 + k_12 / 2))
    k_31 = h * (u_2 + k_22 / 2)
    # ! CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_32 = h * (4 - (1 / 4) * pow(x + h / 2, 3) - (u_1 + k_21 / 2) * (u_2 + k_22 / 2))
    k_41 = h * (u_2 + k_32)
    # ! CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_42 = h * (4 - (1 / 4) * pow(x + h, 3) - (u_1 + k_31) * (u_2 + k_32))

    # ? New values for u_1 and u_2
    new_u_1 = u_1 + (1.0 / 6.0) * (k_11 + 2 * k_21 + 2 * k_31 + k_41)
    new_u_2 = u_2 + (1.0 / 6.0) * (k_12 + 2 * k_22 + 2 * k_32 + k_42)

    row = np.round(
        np.array(
            [
                (u_1, u_2),
                (u_2, u_2_eval),
                (k_11, k_12),
                (k_21, k_22),
                (k_31, k_32),
                (k_41, k_42)
            ]),
        6)

    return new_u_1, new_u_2, row


def generate_kutta_table(start_iteration, iterations, step, u_1, u_2):
    table = []
    for itr in range(0, iterations):
        new_values = runge_kutta_fourth_order(
            x=start_iteration,
            u_1=u_1,
            u_2=u_2,
            h=step
        )
        start_iteration += step
        u_1 = new_values[0]
        u_2 = new_values[1]
        table.append(new_values[2])

    return table


def interpolate(table, aprox_value):
    return round(
        ((aprox_value - table[-1][0]) / (table[-2][0] - table[-1][0])) * table[-2][1] +
        ((aprox_value - table[-2][0]) / (table[-1][0] - table[-2][0])) * table[-1][1]
        , 6)


def linear_shooting(start_iteration, iterations, step, u_1, u_2_first, u_2_second, aprox_value, tolerance):
    runge_kutta_tables = []
    interpolation_table = []

    # * Initial steps (generating two tables to compute an interpolation)
    first_try = generate_kutta_table(start_iteration, iterations, step, u_1, u_2_first)
    second_try = generate_kutta_table(start_iteration, iterations, step, u_1, u_2_second)
    runge_kutta_tables.append(first_try)
    runge_kutta_tables.append(second_try)
    interpolation_table.append([first_try[-1][0][0], u_2_first])
    interpolation_table.append([second_try[-1][0][0], u_2_second])

    # * Interpolation calculation and kutta table
    interpolation = interpolate(interpolation_table, aprox_value)
    next_try = generate_kutta_table(start_iteration, iterations, step, u_1, interpolation)
    runge_kutta_tables.append(next_try)
    interpolation_table.append([next_try[-1][0][0], interpolation])

    # * Compute another interpolation value until difference is less than the tolerance or more than 10 iterations
    itr_counter = 0
    while abs(next_try[-1][0][0] - aprox_value) >= tolerance or itr_counter >= 10:
        interpolation = interpolate(interpolation_table, aprox_value)
        next_try = generate_kutta_table(start_iteration, iterations, step, u_1, interpolation)
        runge_kutta_tables.append(next_try)
        interpolation_table.append([next_try[-1][0][0], interpolation])
        itr_counter += 1

    return next_try[-1][0][0], runge_kutta_tables, interpolation_table


if __name__ == "__main__":
    print("Metodo de disparo lineal")
    u_1 = float(input("Ingresa el valor de y(x): "))
    start_iteration = int(input("Ingresa el valor inicial x: "))
    u_2_first = float(input("Ingresa una aproximacion para y'(x): "))
    u_2_second = float(input("Ingresa otra aproximacion para y'(x): "))
    last_iteration = int(input("Ingresa el valor final x: "))
    aprox_value = float(input("Ingresa el valor a aproximar: "))
    tolerance = float(input("Ingresa el valor de la tolerancia: "))
    step = float(input("Ingresa el valor de h: "))

    headers = np.array(["u/u'", "f(u)/f(u')", "k1", "k2", "k3", "k4"])
    headers_interpolation = np.array(["x", "x'"])

    # ? Adding one extra to include initial iteration
    iterations = int(abs(last_iteration - start_iteration) / step) + 1
    inter, runge_kutta_tables, interpolation_table = linear_shooting(start_iteration, iterations,
                                                                     step, u_1, u_2_first, u_2_second,
                                                                     aprox_value, tolerance)
    print("TABLAS DE RUNGE KUTTA")
    for table in runge_kutta_tables:
        print(tabulate(tabular_data=table, headers=headers) + "\n")

    print("TABLA DE INTERPOLACION")
    print(tabulate(tabular_data=interpolation_table, headers=headers_interpolation))
