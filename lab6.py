from math import sqrt, sin, cos, atan, tan, pi, degrees, radians, asin, atan2, floor
from typing import Tuple, NewType
from functools import wraps

# выводит информация о затраченной оперативной памяти
from memory_profiler import profile
# profile - декортатор!

# Исходные данные
a = 6_378_245
alpha = 1 / 298.3
b = a * (1 - alpha)
e2 = alpha * (2 - alpha)
ee2 = e2 / (1 - e2)
B1 = radians(41 + 11 / 60 + 48.32 / 3600)
A1 = radians(55)
L1 = radians(21)
radian = NewType('radian', float)
degree = NewType('degree', float)
S = 275000
print('Начальные:')
print('Широта b1: ', degrees(B1))
print('Долгота l1: ', degrees(L1))
print('Азимут a1: ', degrees(A1))


def print_dms(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if isinstance(res, tuple):
            print(func.__name__)
            for r in res:
                degree = floor(r)
                mins = (r - degree) * 60
                sec = (mins - floor(mins)) * 60
                print(f"degree: {degree}, mins: {floor(mins)}, sec: {sec:.12f}")
            print('\n')
        return res
    return wrapper


def print_error(*args):
    print(f"dB1={(args[0] - args[1]) * 3600:.10f},  dL1={(args[2] - args[3]) * 3600:.10f}, dA1={(args[4] - abs(args[5] - 180)) * 3600:.10f}")


@print_dms
def direct_task(B1: radian, L1: radian, A1: radian, S: int) -> Tuple[degree, degree, degree]:
    W1 = sqrt(1 - e2 * sin(B1) ** 2)
    sin_u = sin(B1) * sqrt(1 - e2) / W1
    cos_u = cos(B1) / W1
    sin_A0 = cos_u * sin(A1)
    cat_sigma = cos_u * cos(A1) / sin_u
    sin_2sig = 2 * cat_sigma / (cat_sigma ** 2 + 1)
    cos_2sig = (cat_sigma ** 2 - 1) / (cat_sigma ** 2 + 1)
    A = 6_356_863.020 + (10_708.949 - 13.474 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    B = (5354.469 - 8.978 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    C = (2.238 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2) + 0.006
    alpha = 691.46768 - (0.58143 - 0.00144 *
                         (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    betta = (0.2907 - 0.0010 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    sigma0 = (S - (B + C * cos_2sig) * sin_2sig) / A
    sin_2ss0 = sin_2sig * cos(2 * sigma0) + cos_2sig * sin(2 * sigma0)
    cos_2ss0 = cos_2sig * cos(2 * sigma0) - sin_2sig * sin(2 * sigma0)
    sigma = sigma0 + (B + 5 * C * cos_2ss0) * sin_2ss0 / A
    delta = (alpha * sigma + betta * (sin_2ss0 - sin_2sig)) * sin_A0
    delta = radians(delta / 3600)
    sin_u2 = sin_u * cos(sigma) + cos_u * cos(A1) * sin(sigma)
    B2 = atan(sin_u2 / (sqrt(1 - e2) * sqrt(1 - sin_u2 ** 2)))
    lam = atan(sin(A1) * sin(sigma) /
               (cos_u * cos(sigma) - sin_u * sin(sigma) * cos(A1)))
    if sin(A1) >= 0 and tan(lam) >= 0:
        lam = abs(lam)
    elif sin(A1) >= 0 and tan(lam) < 0:
        lam = pi - abs(lam)
    elif sin(A1) < 0 and tan(lam) < 0:
        lam = -abs(lam)
    elif sin(A1) < 0 and tan(lam) >= 0:
        lam = abs(lam) - pi
    L2 = L1 + lam - delta
    A2 = atan(cos_u * sin(A1) / (cos_u * cos(sigma)
                                 * cos(A1) - sin_u * sin(sigma)))
    if sin(A1) < 0 and tan(lam) >= 0:
        A2 = abs(A2)
    elif sin(A1) < 0 and tan(lam) < 0:
        A2 = pi - abs(A2)
    elif sin(A1) >= 0 and tan(lam) >= 0:
        A2 = pi + abs(A2)
    elif sin(A1) >= 0 and tan(lam) > 0:
        A2 = 2 * pi - abs(A2)
    return degree(degrees(B2)), degree(degrees(L2)), degree(degrees(A2))


@print_dms
def reverse_task(B1: radian, L1: radian, B2: radian, L2: radian) -> Tuple[int, degree]:
    l = L2 - L1  # noqa: E741
    cos_B1 = cos(B1)
    sin_B1 = sin(B1)
    cos_B2 = cos(B2)
    sin_B2 = sin(B2)
    sin_l = sin(l)
    cos_l = cos(l)
    p = cos_B2 * sin_l
    q = cos_B1 * sin_B2 - sin_B1 * cos_B2 * cos_l
    tan_a1 = p / q
    a1 = atan2(p, q)
    sin_s = p * sin(a1) + q * cos(a1)
    sigma = asin(sin_s)
    cos_s = sin_B1 * sin_B2 + cos_B1 * cos_B2 * cos_l
    sin_a0 = cos_B1 * sin(a1)
    cos_2_a0 = 1 - sin_a0 ** 2
    A = 6_378_245.0 - (10_601.6 + 84.8 * cos_2_a0) * cos_2_a0
    Bb = 31_947.8 + 125.0 * cos_2_a0
    Cc = 67.0
    x = 2 * sin_B1 * sin_B2 - cos_2_a0 * cos_s
    y = (cos_2_a0 ** 2 - 2 * x ** 2) * cos_s
    S = A * sigma + (Bb * x - Cc * y) * sin_s
    A1 = atan((1 + ee2 * cos_B1 ** 2) * tan_a1)
    return S, degree(degrees(A1))


def funcs(B, L, A):
    M = a * (1 - e2) / pow((1 - e2*sin(B)**2), 3/2)
    N = a * pow((1 - e2 * sin(B)**2), -1/2)
    f_b = cos(A) / M
    f_l = sin(A) / (N * cos(B))
    f_a = sin(A) * tan(B) / N
    return f_b, f_l, f_a


@print_dms
def runge_kutta(h: int, B1: radian, L1: radian, A1: radian, S: int) -> Tuple[degree, degree, degree]:
    nstep = S / h
    B0 = B1
    L0 = L1
    A0 = A1
    for i in range(int(nstep)):
        k1_b = funcs(B0, L0, A0)[0]
        k1_l = funcs(B0, L0, A0)[1]
        k1_a = funcs(B0, L0, A0)[2]
        k2_b = funcs(B0 + h / 2 * k1_b, L0 + h /
                     2 * k1_l, A0 + h / 2 * k1_a)[0]
        k2_l = funcs(B0 + h / 2 * k1_b, L0 + h /
                     2 * k1_l, A0 + h / 2 * k1_a)[1]
        k2_a = funcs(B0 + h / 2 * k1_b, L0 + h /
                     2 * k1_l, A0 + h / 2 * k1_a)[2]
        k3_b = funcs(B0 + h / 2 * k2_b, L0 + h /
                     2 * k2_l, A0 + h / 2 * k2_a)[0]
        k3_l = funcs(B0 + h / 2 * k2_b, L0 + h /
                     2 * k2_l, A0 + h / 2 * k2_a)[1]
        k3_a = funcs(B0 + h / 2 * k2_b, L0 + h /
                     2 * k2_l, A0 + h / 2 * k2_a)[2]
        k4_b = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[0]
        k4_l = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[1]
        k4_a = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[2]
        B0 = B0 + h / 6 * (k1_b + 2 * k2_b + 2 * k3_b + k4_b)
        L0 = L0 + h / 6 * (k1_l + 2 * k2_l + 2 * k3_l + k4_l)
        A0 = A0 + h / 6 * (k1_a + 2 * k2_a + 2 * k3_a + k4_a)
    return degree(degrees(B0)), degree(degrees(L0)), degree(degrees(A0))


@print_dms
def runge_merson(B1: radian, L1: radian, A1: radian, S: int)-> Tuple[degree, degree, degree]:
    eps = 1e-11
    h = 1
    B0 = B1
    L0 = L1
    A0 = A1
    step = 0
    while step <= S:
        k1_b = h / 3 * funcs(B0, L0, A0)[0]
        k1_l = h / 3 * funcs(B0, L0, A0)[1]
        k1_a = h / 3 * funcs(B0, L0, A0)[2]
        k2_b = h / 3 * funcs(B0 + k1_b, L0 + k1_l, A0 + k1_a)[0]
        k2_l = h / 3 * funcs(B0 + k1_b, L0 + k1_l, A0 + k1_a)[1]
        k2_a = h / 3 * funcs(B0 + k1_b, L0 + k1_l, A0 + k1_a)[2]
        k3_b = h / 3 * funcs(
            B0 + 1/2 * k1_b + 1/2 * k2_b,
            L0 + 1/2 * k1_l + 1/2 * k2_l,
            A0 + 1/2 * k1_a + 1/2 * k2_a
        )[0]
        k3_l = h / 3 * funcs(
            B0 + 1/2 * k1_b + 1/2 * k2_b,
            L0 + 1/2 * k1_l + 1/2 * k2_l,
            A0 + 1/2 * k1_a + 1/2 * k2_a
        )[1]
        k3_a = h / 3 * funcs(
            B0 + 1/2 * k1_b + 1/2 * k2_b,
            L0 + 1/2 * k1_l + 1/2 * k2_l,
            A0 + 1/2 * k1_a + 1/2 * k2_a
        )[2]
        k4_b = h / 3 * funcs(
            B0 + 3/8 * k1_b + 9/8 * k3_b,
            L0 + 3/8 * k1_l + 9/8 * k3_l,
            A0 + 3/8 * k1_a + 9/8 * k3_a
        )[0]
        k4_l = h / 3 * funcs(
            B0 + 3/8 * k1_b + 9/8 * k3_b,
            L0 + 3/8 * k1_l + 9/8 * k3_l,
            A0 + 3/8 * k1_a + 9/8 * k3_a
        )[1]
        k4_a = h / 3 * funcs(
            B0 + 3/8 * k1_b + 9/8 * k3_b,
            L0 + 3/8 * k1_l + 9/8 * k3_l,
            A0 + 3/8 * k1_a + 9/8 * k3_a
        )[2]
        k5_b = h / 3 * funcs(
            B0 + 3/2 * k1_b - 9/2 * k3_b + 6 * k4_b,
            L0 + 3/2 * k1_l - 9/2 * k3_l + 6 * k4_l,
            A0 + 3/2 * k1_a - 9/2 * k3_a + 6 * k4_a
        )[0]
        k5_l = h / 3 * funcs(
            B0 + 3/2 * k1_b - 9/2 * k3_b + 6 * k4_b,
            L0 + 3/2 * k1_l - 9/2 * k3_l + 6 * k4_l,
            A0 + 3/2 * k1_a - 9/2 * k3_a + 6 * k4_a
        )[1]
        k5_a = h / 3 * funcs(
            B0 + 3/2 * k1_b - 9/2 * k3_b + 6 * k4_b,
            L0 + 3/2 * k1_l - 9/2 * k3_l + 6 * k4_l,
            A0 + 3/2 * k1_a - 9/2 * k3_a + 6 * k4_a
        )[2]
        eps_b = 0.2 * k1_b - 0.9 * k3_b + 0.8 * k4_b - 0.1 * k5_b
        eps_l = 0.2 * k1_l - 0.9 * k3_l + 0.8 * k4_l - 0.1 * k5_l
        eps_a = 0.2 * k1_a - 0.9 * k3_a + 0.8 * k4_a - 0.1 * k5_a
        if abs(eps_b) > eps and abs(eps_l) > eps and abs(eps_a) > eps:
            h /= 2
        elif 32 * (eps_b) <= eps and 32 * (eps_l) <= eps and 32 * (eps_a) <= eps:
            h *= 2
        B0 = B0 + 1/2 * (k1_b + 4 * k4_b + k5_b)
        L0 = L0 + 1/2 * (k1_l + 4 * k4_l + k5_l)
        A0 = A0 + 1/2 * (k1_a + 4 * k4_a + k5_a)
        if step == S:
            return degrees(B0), degrees(L0), degrees(A0)
        if step + h > S:
            h = S - step
        step += h


@profile
def main():
    print('\nПрямая задача')
    B2_1, L2_1, A2_1 = direct_task(B1, L1, A1, S)
    print('Вычисление ошибки для прямой задачи: ')
    B1_11, L1_11, A1_11 = direct_task(radians(B2_1), radians(L2_1), radians(A2_1), S)
    print_error(degrees(B1), B1_11, degrees(L1), L1_11, degrees(A1), A1_11)

    print('\nВычисление S, A1. Обратная задача\n')
    SS, AA1 = reverse_task(B1, L1, radians(B2_1), radians(L2_1))
    print('Обратная, S:', SS)
    print(f'Ошибка прямая: S: {abs(S - SS)},   A1:  {abs(AA1-degrees(A1)) * 3600}  ')

    print('\nКлассический метод Рунге-Кутты: ')
    print('h = ', S)
    B2, L2, A2 = runge_kutta(S, B1, L1, A1, S)
    print('\n')
    print('h = ', 500)
    B2, L2, A2 = runge_kutta(500, B1, L1, A1, S)
    print('\nh = ', 1)
    B2, L2, A2 = runge_kutta(1, B1, L1, A1, S)
    B2_r, L2_r, A2_r = B2, L2, A2
    B1_1, L1_1, A1_1 = runge_kutta(S, radians(B2), radians(L2), radians(180+A2), S)
    print('Ошибка рунге, h = ', 1)
    print_error(degrees(B1), B1_1, degrees(L1), L1_1, degrees(A1), A1_1)

    print('Ошибка ПРЯМОЙ-РУНГЕ:')
    print('Широта: ', (B2_1 - B2)*3600)
    print('Долгота: ', (L2_1 - L2)*3600)
    print('Азимут: ', (abs(A2_1-180) - A2)*3600)

    print('\nОбратная задача, Рунге-Кутты: ')
    SS, AA1 = reverse_task(B1, L1, radians(B2), radians(L2))
    print('Обратная, S:', SS)
    print(f'Ошибка рунге: S: {abs(S - SS)},   A1:  {abs(AA1-degrees(A1)) * 3600}  ')

    print('\n Метод Рунге-Мерсона')
    B2, L2, A2 = runge_merson(B1, L1, A1, S)
    B1_1, L1_1, A1_1 = runge_merson(radians(B2), radians(L2), radians(A2+180), S)
    print('Ошибка мерсона: ')
    print_error(degrees(B1), B1_1, degrees(L1), L1_1, degrees(A1), A1_1)
    
    print('\nОшибка ПРЯМОЙ-МЕРСОНА:')
    print('Широта: ', (B2_1 - B2)*3600)
    print('Долгота: ', (L2_1 - L2)*3600)
    print('Азимут: ', (abs(A2_1-180) - A2)*3600)
    print('\nОшибка РУНГЕ-МЕРСОНА:')
    print('Широта: ', (B2_r - B2)*3600)
    print('Долгота: ', (L2_r - L2)*3600)
    print('Азимут: ', (A2_r - A2)*3600)


if __name__ == '__main__':
    main()
