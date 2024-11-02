for i in range(1000, 10000):
    if is_prime(i):
        for j, _ in factor(i - 1):
            if j > 3:
                break
        else:
            print(i)
