
from scipy.stats import chisquare, norm
from scipy.stats import chi2
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
from collections import Counter
import pandas as pd
import numpy as np
# Generatory

def xorshift128(state, n):
    x, y, z, w = [np.uint32(s) for s in state]
    vector = np.zeros(n)
    for i in range(n):
        t = x ^ (x << np.uint32(11))
        x, y, z = y, z, w
        w = w ^ (w >> np.uint32(19)) ^ t ^ (t >> np.uint32(8))
        w = np.uint32(w)
        vector[i] = int(w) / 2**32
    return vector

def rc4_prng(K, n, m=256):
    L = len(K)
    S = [i for i in range(m)]
    j = 0
    for i in range(m):
        j = (j + S[i] + K[i % L]) % m
        S[i], S[j] = S[j], S[i]

    i = 0
    j = 0
    vector = np.zeros(n)
    for idx in range(n):
        i = (i + 1) % m
        j = (j + S[i]) % m
        S[i], S[j] = S[j], S[i]
        r = S[(S[i] + S[j]) % m]
        vector[idx] = r / m
    return vector

def lcg_to_01(state, n, a=1664525, c=1013904223, m=2**32):
  vector = np.zeros(n)
  for i in range(n):
    state = (a * state + c) % m
    vector[i] = state / m
  return vector

def lcg(state, n, a=1664525, c=1013904223, m=2**32):
  vector = np.zeros(n)
  for i in range(n):
    state = (a * state + c) % m
    vector[i] = state
  return vector

# Histogram wartości dla n=2^12
n = 2**12
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])

k = 20
L = [1, 2, 3, 4, 5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1,1,1], n)



generators = {
    "LCG": vector_LCG,
    "Xorshift": vector_xor,
    "RC4": vector_rc4,
    "Mersenne Twister": vector_MT
}

k =32

plt.figure(figsize=(12, 8))


for i, (name, vector) in enumerate(generators.items(), 1):
    plt.subplot(2, 2, i)
    plt.hist(vector, bins=k, range=(0,1), edgecolor='black')
    plt.title(name)
    plt.xlabel("Bin")
    plt.ylabel("Values count")
    plt.ylim(0, max(len(vector)/k * 1.2, 200))
    plt.axhline(len(vector)/k, color='red', linestyle='dashed', label='Expected number')
    plt.legend()

plt.tight_layout()
plt.show()

# Histogram wartości dla n=50000

n = 50000
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])

k = 20
L = [1, 2, 3, 4, 5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1,1,1], n)

import matplotlib.pyplot as plt


generators = {
    "LCG": vector_LCG,
    "Xorshift": vector_xor,
    "RC4": vector_rc4,
    "Mersenne Twister": vector_MT
}

k =32

plt.figure(figsize=(12, 8))


for i, (name, vector) in enumerate(generators.items(), 1):
    plt.subplot(2, 2, i)
    plt.hist(vector, bins=k, range=(0,1), edgecolor='black')
    plt.title(name)
    plt.xlabel("Bin")
    plt.ylabel("Values count")
    plt.ylim(0, max(len(vector)/k * 1.2, 200))
    plt.axhline(len(vector)/k, color='red', linestyle='dashed', label='Expected number')
    plt.legend()

plt.tight_layout()
plt.show()

# Testy

def runs_test(data):
    threshold = np.median(data)
    Y = np.where(data > threshold, 1, 0)
    n = len(Y)

    R = 1 + np.sum(Y[1:] != Y[:-1])
    n1, n0 = np.sum(Y), n - np.sum(Y)
    E_R = (2 * n1 * n0) / (n + 1)
    Var_R = ((2 * n1 * n0 ) / (n**2 * (n - 1))) * (2 * n1 * n0 - n)
    Z = (R - E_R) / np.sqrt(Var_R)
    p = 2 * (1 - norm.cdf(abs(Z)))
    # print(f"fr={Var_R}")
    # print(f"n= {n}")
    # print((n**2 * (n - 1)))
    return R, E_R, Z, p

def birthday_test(values, m=2**32):
    k = len(values)
    Y = np.array(values, dtype=np.int64)

    Y_sorted = np.sort(Y)
    S = [(Y_sorted[i+1] - Y_sorted[i]) for i in range(k-1)]
    S.append((Y[0] + m) - Y[-1])

    counts = Counter(S)
    counts = np.array(list(counts.values()))
    K = np.sum(counts[counts > 1] - 1)

    lambda1 = k**3 / (4*m)
    p = 1 - poisson.cdf(K - 1, lambda1)

    return K, p, lambda1

def monobit_test(bits):
    bits = np.array(bits)
    n = len(bits)

    x = 2*bits - 1

    s_n = np.sum(x) / np.sqrt(n)

    p_value = 2 * (1 - norm.cdf(abs(s_n)))

    return s_n, p_value


bits = np.zeros(len(vector))

def bits_from_threshold(float_list, threshold=0.5):
    return [1 if x > threshold else 0 for x in float_list]


def chi2_test(vector, k):
    n = len(vector)
    E = n / k
    bins = np.linspace(0, 1, k + 1)

    categories = np.digitize(vector, bins) -1

    Y = np.bincount(categories, minlength=k)

    chi2_stat = np.sum((Y - E)**2 / E)
    df = k - 1
    p_value = 1 - chi2.cdf(chi2_stat, df)

    return chi2_stat, p_value


# Wyniki testu dla n=2^12
n = 2**12
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])

k = 32
L = [1, 2, 3, 4, 5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1,1,1], n)

vector = lcg(1, n)
vector_MT_new = (vector_MT* 2**32).astype(np.uint32)
vector_rc4_new = (vector_rc4* 2**32).astype(np.uint32)
vector_xor_new = (vector_xor* 2**32).astype(np.uint32)



chi2_vals = {
    "LCG": chi2_test(vector_LCG, k)[1],
    "Xorshift128": chi2_test(vector_xor, k)[1],
    "RC4": chi2_test(vector_rc4, k)[1],
    "MT19937": chi2_test(vector_MT, k)[1]
}

mono_vals = {
    "LCG": monobit_test(bits_from_threshold(vector_LCG))[1],
    "Xorshift128": monobit_test(bits_from_threshold(vector_xor))[1],
    "RC4": monobit_test(bits_from_threshold(vector_rc4))[1],
    "MT19937": monobit_test(bits_from_threshold(vector_MT))[1]
}

runs_vals = {
    "LCG": runs_test(vector_LCG)[3],
    "Xorshift128": runs_test(vector_xor)[3],
    "RC4": runs_test(vector_rc4)[3],
    "MT19937": runs_test(vector_MT)[3]
}

vector_MT_new = np.floor(vector_MT * 2**32).astype(np.uint32)
vector_rc4_new = np.floor(vector_rc4 * 2**32).astype(np.uint32)
vector_xor_new = np.floor(vector_xor * 2**32).astype(np.uint32)
vector_LCG_new = np.floor(vector_LCG * 2**32).astype(np.uint32)

birthday_vals = {
    "LCG": birthday_test(vector_LCG_new)[1],
    "Xorshift128": birthday_test(vector_xor_new)[1],
    "RC4": birthday_test(vector_rc4_new)[1],
    "MT19937": birthday_test(vector_MT_new)[1]
}


df = pd.DataFrame({
    "Generator": ["LCG", "Xorshift128", "RC4", "MT19937"],
    "Chi-square": [chi2_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Runs": [runs_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Monobit": [mono_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Birthday Spacing": [birthday_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]]
})

print(df)

# Wyniki testu dla n=50000

n = 50000
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])

k = 32
L = [1, 2, 3, 4, 5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1,1,1], n)

vector = lcg(1, n)
vector_MT_new = (vector_MT* 2**32).astype(np.uint32)
vector_rc4_new = (vector_rc4* 2**32).astype(np.uint32)
vector_xor_new = (vector_xor* 2**32).astype(np.uint32)


chi2_vals = {
    "LCG": chi2_test(vector_LCG, k)[1],
    "Xorshift128": chi2_test(vector_xor, k)[1],
    "RC4": chi2_test(vector_rc4, k)[1],
    "MT19937": chi2_test(vector_MT, k)[1]
}

mono_vals = {
    "LCG": monobit_test(bits_from_threshold(vector_LCG))[1],
    "Xorshift128": monobit_test(bits_from_threshold(vector_xor))[1],
    "RC4": monobit_test(bits_from_threshold(vector_rc4))[1],
    "MT19937": monobit_test(bits_from_threshold(vector_MT))[1]
}

runs_vals = {
    "LCG": runs_test(vector_LCG)[3],
    "Xorshift128": runs_test(vector_xor)[3],
    "RC4": runs_test(vector_rc4)[3],
    "MT19937": runs_test(vector_MT)[3]
}

vector_MT_new = np.floor(vector_MT * 2**32).astype(np.uint32)
vector_rc4_new = np.floor(vector_rc4 * 2**32).astype(np.uint32)
vector_xor_new = np.floor(vector_xor * 2**32).astype(np.uint32)
vector_LCG_new = np.floor(vector_LCG * 2**32).astype(np.uint32)

birthday_vals = {
    "LCG": birthday_test(vector_LCG_new)[1],
    "Xorshift128": birthday_test(vector_xor_new)[1],
    "RC4": birthday_test(vector_rc4_new)[1],
    "MT19937": birthday_test(vector_MT_new)[1]
}


df = pd.DataFrame({
    "Generator": ["LCG", "Xorshift128", "RC4", "MT19937"],
    "Chi-square": [chi2_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Runs": [runs_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Monobit": [mono_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]],
    "Birthday Spacing": [birthday_vals[g] for g in ["LCG", "Xorshift128", "RC4", "MT19937"]]
})

print(df)

# heatmapa dla n=2^12

n = 2**12
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])
L = [1,2,3,4,5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1, 1, 1], n)

generators = {
    "LCG": vector_LCG,
    "Xorshift128": vector_xor,
    "RC4": vector_rc4,
    "MT19937": vector_MT
}



fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for ax, (name, data) in zip(axes, generators.items()):
    h = ax.hist2d(
        data[:-1], data[1:], bins=30,
        range=[[0, 1], [0, 1]],
        cmap='viridis'
    )
    ax.set_title(name)
    ax.set_xlabel('$x_i$')
    ax.set_ylabel('$x_{i+1}$')
    fig.colorbar(h[3], ax=ax, label="Values count")

plt.tight_layout()
plt.show()

# heatmapa dla n=50000

n = 50000
random.seed(5678)
vector_MT = np.array([random.random() for _ in range(n)])
L = [1,2,3,4,5]
vector_LCG = lcg_to_01(1, n)
vector_rc4 = rc4_prng(L, n)
vector_xor = xorshift128([1, 1, 1, 1], n)

generators = {
    "LCG": vector_LCG,
    "Xorshift128": vector_xor,
    "RC4": vector_rc4,
    "MT19937": vector_MT
}



fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for ax, (name, data) in zip(axes, generators.items()):
    h = ax.hist2d(
        data[:-1], data[1:], bins=30,
        range=[[0, 1], [0, 1]],
        cmap='viridis'
    )
    ax.set_title(name)
    ax.set_xlabel('$x_i$')
    ax.set_ylabel('$x_{i+1}$')
    fig.colorbar(h[3], ax=ax, label="Values count")

plt.tight_layout()
plt.show()

# Second-level testing

def second_level_test(vector_random, num_blocks=500, k_bins=10, k_second=10):
    n = len(vector_random)
    block_size = n // num_blocks
    p_values = []

    for i in range(num_blocks):
        block = vector_random[i * block_size:(i + 1) * block_size]
        _, p = chi2_test(block, k_bins)
        p_values.append(p)

    p_values = np.array(p_values)
    counts, _ = np.histogram(p_values, bins=np.linspace(0, 1, k_second + 1))
    expected = len(p_values) / k_second
    chi2_stat_2 = np.sum((counts - expected)**2 / expected)
    p_second = 1 - chi2.cdf(chi2_stat_2, df=k_second - 1)

    return p_values, chi2_stat_2, p_second

n=50000

def compare_generators_plots(n=50000, num_blocks=500, k_bins=10, k_second=10):
    random.seed(5678)
    vector_MT = np.array([random.random() for _ in range(n)])
    vector_XOR = xorshift128([1,1,1,1], n)
    vector_RC4 = rc4_prng([1, 2, 3, 4, 5], n)
    vector_LCG = lcg_to_01(1, n)

    gens = {
        "LCG": vector_LCG,
        "Xorshift128": vector_XOR,
        "RC4": vector_RC4,
        "MT19937": vector_MT
    }

    results = {}
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    for ax, (name, vec) in zip(axes, gens.items()):
        p_values, chi2_stat_2, p_second = second_level_test(vec, num_blocks, k_bins, k_second)
        mean_p = np.mean(p_values)
        results[name] = (p_second, mean_p)

        counts, _ = np.histogram(p_values, bins=np.linspace(0, 1, k_second + 1))
        expected = len(p_values) / k_second

        ax.bar(np.linspace(0.05, 0.95, k_second), counts, width=0.09, edgecolor='black')
        ax.axhline(expected, color='red', linestyle='--', label='Expected value')
        ax.set_title(f"{name}")
        ax.set_ylim(0, max(counts) * 1.2)
        ax.set_xlabel("p-value")
        ax.set_ylabel("Values count")
        ax.legend()

    plt.tight_layout()
    plt.show()

    print("\n=== PODSUMOWANIE (p-value 2nd level) ===")
    for gen, (p2, meanp) in results.items():
        print(f"{gen:10s}: p(2nd) = {p2:.3f},  mean p(1st) = {meanp:.3f}")

    return results


compare_generators_plots()