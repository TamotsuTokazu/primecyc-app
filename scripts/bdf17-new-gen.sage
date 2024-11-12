p0, p1 = 73, 97
assert is_prime(p0) and is_prime(p1)

n = max(nn for nn in [1, 10, 100, 600] if nn < p0 - 2)

t = 2 ** 6

m = lcm(p0 * p1 * (p0 - 1) * (p1 - 1), t)

N = 2 ** 59
N = (N // m) * m + 1
while not is_prime(N):
    N += m

print(
f'''const Integer t("{t}");
const usint n = {n};
const usint p0 = {p0};
const usint p1 = {p1};
const usint pq = p0 * p1;
const usint Bks = 1 << 8;
const Integer Q("{N}");
''')

g = primitive_root(N)

print(
f'''const Integer rootOfUnity("{pow(g, (N - 1) // (p0 * (p0 - 1) * p1 * (p1 - 1)), N)}");
const Integer rootOfUnity0("{pow(g, (N - 1) // (p0 * (p0 - 1)), N)}");
const Integer rootOfUnity1("{pow(g, (N - 1) // (p1 * (p1 - 1)), N)}");
''')