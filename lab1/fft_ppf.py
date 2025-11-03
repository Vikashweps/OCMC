import numpy as np
import time
import matplotlib.pyplot as plt

f = 5
T = 1/f
Ts = 0.001
Ns = [128, 256, 512, 1024, 2048, 4096, 2**16]
iters = 200

t_ppf = []
t_fft = []

for N in Ns:
    t = np.arange(0, N * Ts, Ts)
    s = 2 * np.cos(2 * np.pi * f * t)

    sum_ppf = 0.0
    sum_fft = 0.0

    for _ in range(iters):
        # Прямой метод
        start = time.perf_counter()
        for n in range(5):
            c = np.cos(2 * np.pi * n * f * t)
            sc = np.sin(2 * np.pi * n * f * t)
            if n == 0:
                1/T*np.trapezoid(s * c, t)
                1/T*np.trapezoid(s * sc, t)
            else:
                2/T*np.trapezoid(s * c, t)
                2/T*np.trapezoid(s * sc, t)
        sum_ppf += time.perf_counter() - start

        # FFT
        start = time.perf_counter()
        F = np.fft.fft(s)

        sum_fft += time.perf_counter() - start

    t_ppf.append(sum_ppf / iters)
    t_fft.append(sum_fft / iters)

plt.figure(figsize=(8, 5))
plt.semilogy(Ns, t_ppf, 'o-', label='Прямой метод')
plt.semilogy(Ns, t_fft, 'o-', label='FFT')
plt.xlabel('N')
plt.ylabel('Время (сек), лог. шкала')
plt.title(f'Производительность ')
plt.legend()
plt.grid(True)
plt.xticks(Ns)
plt.tight_layout()
plt.show()