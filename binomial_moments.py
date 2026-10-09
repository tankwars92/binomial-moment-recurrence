from math import comb


def alternating_binomial_moment(n, m):
    if n < 0 or m < 0:
        raise ValueError("n and m must be non-negative")

    if n % 2 == 0:
        f0 = (-1) ** (n // 2) * comb(n, n // 2)
        f1 = (n // 2) * f0
    else:
        q = (n + 1) // 2
        f0 = 0
        f1 = (-1) ** q * q * comb(n, (n - 1) // 2)

    f = [f0]
    if m >= 1:
        f.append(f1)

    for r in range(m - 1):
        numerator = ((3 * r + 2 - 2 * n) * f[r + 1]
                     + (n - r) ** 2 * f[r])
        if numerator % 2:
            raise ArithmeticError("Recurrence divisibility failure")
        f.append(-numerator // 2)

    stirling = [0] * (m + 1)
    stirling[0] = 1
    for _ in range(m):
        row = [0] * (m + 1)
        for j in range(1, m + 1):
            row[j] = stirling[j - 1] + j * stirling[j]
        stirling = row

    return sum(stirling[j] * f[j] for j in range(m + 1))


if __name__ == "__main__":
    for n, m in ((0, 0), (5, 3), (12, 2), (100, 10)):
        print("n=%d m=%d A(n,m)=%d" %
              (n, m, alternating_binomial_moment(n, m)))
