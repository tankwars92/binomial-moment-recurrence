from math import comb
from statistics import median
from time import perf_counter
from binomial_moments import alternating_binomial_moment


def direct_sum(n, m):
    total = 0
    c = 1
    for k in range(n + 1):
        term = c * c * (k ** m)
        total += -term if k % 2 else term
        if k < n:
            c = c * (n - k) // (k + 1)
    return total


def krawtchouk_center(n, j):
    total = 0
    for a in range(n % 2, 2 * j + 1, 2):
        b = (n - a) // 2
        if 0 <= b <= n - j:
            term = comb(2 * j, a) * comb(n - j, b)
            total += -term if b % 2 else term
    return total


def krawtchouk_method(n, m):
    polys = [[1]]
    if m >= 1:
        polys.append([n, -2])
    for j in range(1, m):
        p = polys[j]
        prev = polys[j - 1]
        q = [0] * (j + 2)
        for i, c in enumerate(p):
            q[i] += n * c
            q[i + 1] -= 2 * c
        for i, c in enumerate(prev):
            q[i] -= j * (n - j + 1) * c
        polys.append(q)

    moments = []
    for j in range(m + 1):
        rhs = comb(n, j) * krawtchouk_center(n, j)
        rhs *= factorial_int(j)
        lower = sum(polys[j][r] * moments[r] for r in range(j))
        numerator = rhs - lower
        leading = polys[j][j]
        if numerator % leading:
            raise ArithmeticError("Krawtchouk recurrence divisibility failure")
        moments.append(numerator // leading)
    return moments[m]


def factorial_int(n):
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result


def measure(fn, n, m, repeats=5):
    values, timings = [], []
    for _ in range(repeats):
        start = perf_counter()
        values.append(fn(n, m))
        timings.append(perf_counter() - start)
    if any(v != values[0] for v in values):
        raise ArithmeticError("The same method returned different results")
    return values[0], median(timings), timings


def main():
    cases = ((100, 10), (1000, 20), (10000, 20), (10000, 50))
    repeats = 5
    methods = (
        ("Our recurrence", alternating_binomial_moment),
        ("Krawtchouk method", krawtchouk_method),
        ("Direct sum", direct_sum),
    )
    print("Exact comparison of three methods")
    print("A(n,m) = sum((-1)^k * C(n,k)^2 * k^m)")
    print("Each method: %d timed runs; median is reported." % repeats)
    print("The direct sum updates C(n,k) incrementally.")
    for n, m in cases:
        results = {}
        print("\nn=%d, m=%d" % (n, m))
        for label, fn in methods:
            value, elapsed, _ = measure(fn, n, m, repeats)
            results[label] = (value, elapsed)
            print("%-20s %.6f s" % (label + ":", elapsed))
        values = [item[0] for item in results.values()]
        if not all(value == values[0] for value in values):
            raise ArithmeticError("Methods disagree for n=%d, m=%d" % (n, m))
        our_t = results["Our recurrence"][1]
        kraw_t = results["Krawtchouk method"][1]
        direct_t = results["Direct sum"][1]
        print("Exact results match: yes")
        print("Direct / our recurrence: %.2f times" % (direct_t / our_t))
        print("Krawtchouk / our recurrence: %.2f times" % (kraw_t / our_t))
    for n in range(31):
        for m in range(16):
            expected = direct_sum(n, m)
            if alternating_binomial_moment(n, m) != expected:
                raise ArithmeticError("Our recurrence failed small-input check")
            if krawtchouk_method(n, m) != expected:
                raise ArithmeticError("Krawtchouk method failed small-input check")
    print("\nExhaustive check passed for 0 <= n <= 30 and 0 <= m <= 15.")


if __name__ == "__main__":
    main()
