from tabulate import tabulate
import warnings
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)


def u_2_equation(x, u_1, u_2):
    # TODO: CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    return 2 * pow(u_1, 3)


def runge_kutta_fourth_order(x, u_1, u_2, h, table):
    u_1_eval = u_2
    u_2_eval = u_2_equation(x, u_1, u_2)

    k_11 = h * u_1_eval
    k_12 = h * u_2_eval
    k_21 = h * (u_2 + k_12 / 2)
    # TODO: CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_22 = h * 2 * pow((u_1 + k_11 / 2), 3)
    k_31 = h * (u_2 + k_22 / 2)
    # TODO: CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_32 = h * 2 * pow((u_1 + k_21 / 2), 3)
    k_41 = h * (u_2 + k_32)
    # TODO: CHANGE IN CASE FUNCTION DEFINITIONS IS DIFFERENT
    k_42 = h * 2 * pow((u_1 + k_31), 3)

    # New values for u_1 and u_2
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

    table.append(row)

    return new_u_1, new_u_2


if __name__ == "__main__":
    print("Metodo de disparo lineal")
    u_1 = float(input("Ingresa el valor de y(x): "))
    start_iteration = int(input("Ingresa el valor inicial x: "))
    u_2 = float(input("Ingresa una aproximacion para y'(x): "))
    last_iteration = int(input("Ingresa el valor final x: "))
    step = float(input("Ingresa el valor de h: "))

    headers = np.array(["u/u'", "f(u)/f(u')", "k1", "k2", "k3", "k4"])

    runge_kutta_table = []
    # adding one extra to include initial iteration
    iterations = int(abs(last_iteration - start_iteration) / step) + 1
    for itr in range(0, iterations):
        new_values = runge_kutta_fourth_order(
            x=start_iteration,
            u_1=u_1,
            u_2=u_2,
            h=step,
            table=runge_kutta_table
        )
        start_iteration += step
        u_1 = new_values[0]
        u_2 = new_values[1]

    print(tabulate(tabular_data=runge_kutta_table, headers=headers))
