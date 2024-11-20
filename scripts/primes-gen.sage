for i in range(10, 100000):
    if is_prime(i):
        for j, _ in factor(i - 1):
            if j > 3:
                break
        else:
            print(i, factor(i - 1))
