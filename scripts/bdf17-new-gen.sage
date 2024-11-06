p0, p1 = 1153, 1297
assert is_prime(p0) and is_prime(p1)

n = max(nn for nn in [1, 10, 100, 600] if nn < p0 - 2)

t = 2 ** 6 + 1

m = lcm(p0 * p1 * (p0 - 1) * (p1 - 1), t)

N = 2 ** 56
N = (N // m) * m + 1
while not is_prime(N):
    N += m

print(
f'''
const Integer t("{t}");
const uint32_t n = {n};
const uint32_t p0 = {p0};
const uint32_t p1 = {p1};
const uint32_t pq = p0 * p1;
const uint32_t Bks = 1 << 6;
const Integer Q("{N}");
''')

g = primitive_root(N)

print(
f'''
const Integer rootOfUnity("{pow(g, (N - 1) // (p0 * (p0 - 1) * p1 * (p1 - 1)), N)}");
const Integer rootOfUnity0("{pow(g, (N - 1) // (p0 * (p0 - 1)), N)}");
const Integer rootOfUnity1("{pow(g, (N - 1) // (p1 * (p1 - 1)), N)}");
''')