import sys

N = Integer(sys.argv[1])

Qbound = 2 ^ 60
Q = Qbound // N * N + 1
while not is_prime(Q):
    Q -= N

print(Q)
print(pow(primitive_root(Q), (Q - 1) // (N * (N - 1)), Q))