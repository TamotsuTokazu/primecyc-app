p0, p1 = 1153, 1297
assert is_prime(p0) and is_prime(p1)

n = max(nn for nn in [1, 10, 100, 600] if nn < p0 - 2)

t = 2 ** 6 + 1

m = lcm(p0 * p1 * 2, t)

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
const Integer rootOfUnitypq("{pow(g, (N - 1) // (2 * p0 * p1), N)}");
const Integer rootOfUnity0("{pow(g, (N - 1) // (2 * p0), N)}");
const Integer rootOfUnity1("{pow(g, (N - 1) // (2 * p1), N)}");
''')

dftsizep0 = 2 ^ (ceil(log(2 * p0 - 1, 2)))
dftsizep1 = 2 ^ (ceil(log(2 * p1 - 1, 2)))
dftsize = 2 ^ (ceil(log(2 * p0 * p1 - 1, 2)))

N = 2 ^ (10 + ceil(log(2 * p0 * p1 - 1, 2) + 2 * log(N, 2)))
N = (N // dftsize) * dftsize + 1
while not is_prime(N):
    N += dftsize

g = primitive_root(N)

print(
f'''
const Integer bigModulus("{N}");
const Integer bigRootOfUnity0("{pow(g, (N - 1) // dftsizep0, N)}");
const Integer bigRootOfUnity1("{pow(g, (N - 1) // dftsizep1, N)}");
const Integer bigRootOfUnitypq("{pow(g, (N - 1) // dftsize, N)}");
''')
