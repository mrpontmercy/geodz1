from math import sqrt, sin, cos, atan, tan, pi, degrees, radians, asin, atan2, floor
from typing import Tuple, NewType

from memory_profiler import profile

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
S = 35000


def print_dms(func):
    '''Выводит значения float в градусах, минутах, секундах'''
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if isinstance(res, tuple):
            print(func.__name__)
            for r in res:
                degree = floor(r)
                mins = (r - degree) * 60
                sec = (mins - floor(mins)) * 60
                print(f"degree: {degree}, mins: {floor(mins)}, sec: {sec:.3f}")
        return res
    return wrapper


@print_dms
def direct_task(B1: radian, L1: radian, A1: radian, S: radian) -> Tuple[degree, degree, degree]:
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
    alpha = 691.46768 - (0.58143 - 0.00144 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    betta = (0.2907 - 0.0010 * (1 - sin_A0 ** 2)) * (1 - sin_A0 ** 2)
    sigma0 = (S - (B + C * cos_2sig) * sin_2sig) / A
    sin_2ss0 = sin_2sig * cos(2 * sigma0) + cos_2sig * sin(2 * sigma0)
    cos_2ss0 = cos_2sig * cos(2 * sigma0) - sin_2sig * sin(2 * sigma0)
    sigma = sigma0 + (B + 5 * C * cos_2ss0) * sin_2ss0 / A
    delta = (alpha * sigma + betta * (sin_2ss0 - sin_2sig)) * sin_A0
    delta = radians(delta / 3600)
    sin_u2 = sin_u * cos(sigma) + cos_u * cos(A1) * sin(sigma)
    B2 = atan(sin_u2 / (sqrt(1 - e2) * sqrt(1 - sin_u2 ** 2)))
    lam = atan(sin(A1) * sin(sigma) / (cos_u * cos(sigma) - sin_u * sin(sigma) * cos(A1)))
    if sin(A1) >= 0 and tan(lam) >= 0:
        lam = abs(lam)
    elif sin(A1) >= 0 and tan(lam) < 0:
        lam = pi - abs(lam)
    elif sin(A1) < 0 and tan(lam) < 0:
        lam = -abs(lam)
    elif sin(A1) < 0 and tan(lam) >= 0:
        lam = abs(lam) - pi
    L2 = L1 + lam - delta
    A2 = atan(cos_u * sin(A1) / (cos_u * cos(sigma) * cos(A1) - sin_u * sin(sigma)))
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
    cos_B1 = cos(B1); sin_B1 = sin(B1); cos_B2 = cos(B2);
    sin_B2 = sin(B2); sin_l = sin(l); cos_l = cos(l);
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

@print_dms
def runge_kutta(h: int, B1: radian, L1: radian, A1: radian, S: int) -> Tuple[degree, degree, degree]:
    print(h)
    def funcs(B, L, A):
        M = a * (1 - e2) / pow((1 - e2*sin(B)**2), 3/2)
        N = a * pow((1 - e2 * sin(B)**2), -1/2)
        f_b = cos(A) / M
        f_l = sin(A) / (N * cos(B))
        f_a = sin(A) * tan(B) / N
        return f_b, f_l, f_a
    nstep = S / h
    B0 = B1; L0 = L1; A0 = A1
    for i in range(int(nstep)):
        k1_b = funcs(B0, L0, A0)[0]
        k1_l = funcs(B0, L0, A0)[1]
        k1_a = funcs(B0, L0, A0)[2]
        k2_b = funcs(B0 + h / 2 * k1_b, L0 + h / 2 * k1_l, A0 + h / 2 * k1_a)[0]
        k2_l = funcs(B0 + h / 2 * k1_b, L0 + h / 2 * k1_l, A0 + h / 2 * k1_a)[1]
        k2_a = funcs(B0 + h / 2 * k1_b, L0 + h / 2 * k1_l, A0 + h / 2 * k1_a)[2]
        k3_b = funcs(B0 + h / 2 * k2_b, L0 + h / 2 * k2_l, A0 + h / 2 * k2_a)[0]
        k3_l = funcs(B0 + h / 2 * k2_b, L0 + h / 2 * k2_l, A0 + h / 2 * k2_a)[1]
        k3_a = funcs(B0 + h / 2 * k2_b, L0 + h / 2 * k2_l, A0 + h / 2 * k2_a)[2]
        k4_b = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[0]
        k4_l = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[1]
        k4_a = funcs(B0 + h * k3_b, L0 + h * k3_l, A0 + h * k3_a)[2]
        B0 = B0 + h / 6 * (k1_b + 2 * k2_b + 2 * k3_b + k4_b)
        L0 = L0 + h / 6 * (k1_l + 2 * k2_l + 2 * k3_l + k4_l)
        A0 = A0 + h / 6 * (k1_a + 2 * k2_a + 2 * k3_a + k4_a)
    return degree(degrees(B0)), degree(degrees(L0)), degree(degrees(A0))


@profile
def main():
    B2, L2, A2 = direct_task(B1, L1, A1, S)
    SS, AA1 = reverse_task(B1, L1, radians(B2), radians(L2))
    print('Широта b2: ', B2)
    print('Долгота l2: ', L2)
    print('Азимут a2: ', A2)
    print('Обратная задача, S: ', SS)
    print('Обратная задача, A1: ', AA1)
    print('метод рунге-кутты')
    B2, L2, A2 = runge_kutta(S, B1, L1, A1, S)
    print('Широта b2: ', B2)
    print('Долгота l2: ', L2)
    print('Азимут a2: ', A2)


if __name__ == '__main__':
    main()
