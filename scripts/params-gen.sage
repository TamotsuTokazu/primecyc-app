p0 = 1153
p1 = 1297
assert is_prime(p0) and is_prime(p1)

t = 2 ** 6

m = lcm(p0 * p1 * (p0 - 1) * (p1 - 1), t)

N = 2 ** 55
N = (N // m) * m + 1
while not is_prime(N):
    N += m

print(N)

g = primitive_root(N)
print(pow(g, (N - 1) // (p0 * (p0 - 1)), N))
print(pow(g, (N - 1) // (p1 * (p1 - 1)), N))
print(pow(g, (N - 1) // (2 * (p0 - 1) * (p1 - 1)), N))