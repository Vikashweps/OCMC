import numpy as np
import matplotlib.pyplot as plt

# Параметры
sizes = np.linspace(100, 1100, 100)
crc_lengths = [2, 3, 4, 5, 6, 7, 8]
colors = ['darkred', 'red', 'blue', 'green', 'orange', 'purple', 'brown']

plt.figure(figsize=(12, 7))

for crc_idx, crc_len in enumerate(crc_lengths):
    max_prob = 1 / (2 ** crc_len)

    # Логистический рост
    k = 0.008
    N0 = 600
    probs = max_prob / (1 + np.exp(-k * (sizes - N0)))

    plt.plot(sizes, probs, '-',
             color=colors[crc_idx],
             linewidth=2,
             label=f'CRC-{crc_len}')

# Настройки
plt.xlabel('Размер данных (бит)')
plt.ylabel('Вероятность пропуска ошибки')
plt.title('Рост вероятности пропуска ошибок для разных CRC')
plt.legend()
plt.grid(True, alpha=0.3)
plt.ylim(0, 0.3)

plt.tight_layout()
plt.show()