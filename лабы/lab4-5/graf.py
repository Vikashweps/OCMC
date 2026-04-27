import numpy as np
import matplotlib.pyplot as plt

# Функция для вычисления CRC через XOR
def add_crc(data_bits, crc_poly):
    data = data_bits.copy()
    poly_len = len(crc_poly)
    # Добавляем нули для CRC
    data_extended = np.concatenate([data, np.zeros(poly_len - 1, dtype=int)])
    #  XOR
    for i in range(len(data)):
        if data_extended[i] == 1:
            data_extended[i:i + poly_len] ^= crc_poly

    return data_extended[-poly_len + 1:]

np.random.seed(42)
crc_configs = {
    'CRC-2': [1, 1],
    'CRC-3': [1, 0, 1],
    'CRC-4': [1, 0, 0, 1],
    'CRC-5': [1, 0, 1, 0, 1],
    'CRC-6': [1, 1, 0, 0, 0, 1],
    'CRC-7': [1, 0, 0, 1, 0, 0, 1],
    'CRC-8': [1, 1, 0, 1, 1, 0, 1, 0]
}

# Диапазон размеров данных
data_sizes = np.arange(100, 10001, 800, dtype=int)  # от 100 до 10000 с шагом 500
num_trials = 100
colors = ['darkred', 'red', 'blue', 'green', 'orange', 'purple', 'brown']

plt.figure(figsize=(12, 6))

# Для каждого CRC
for (crc_name, crc_poly), color in zip(crc_configs.items(), colors):
    missed_counts = []

    # Для каждого размера данных
    for size in data_sizes:
        missed_errors = 0

        # Испытания
        for _ in range(num_trials):
            # Исходные данные
            original_data = np.random.randint(0, 2, size)

            # Вычисляем CRC
            crc_bits = add_crc(original_data, crc_poly)

            # Передаваемое сообщение
            transmitted = np.concatenate([original_data, crc_bits])

            # Полностью случайное принятое сообщение
            received = np.random.randint(0, 2, len(transmitted))

            # Убедимся, что сообщения не идентичны
            if np.array_equal(received, transmitted):
                error_pos = np.random.randint(0, len(transmitted))
                received[error_pos] ^= 1

            # Проверяем CRC
            received_data = received[:-len(crc_bits)]
            received_crc = received[-len(crc_bits):]
            computed_crc = add_crc(received_data, crc_poly)

            if np.array_equal(received_crc, computed_crc):
                missed_errors += 1

        missed_counts.append(missed_errors)
        # Выводим информацию о текущем размере
        print(f'{crc_name}, размер {size}: пропущено {missed_errors}/{num_trials} ошибок')

    # Построение графика
    plt.plot(data_sizes, missed_counts, 'o-', color=color, linewidth=2,
             markersize=4, label=crc_name)

# Настройки графика
plt.xlabel('Размер данных (бит)')
plt.ylabel('Количество пропущенных ошибок')
plt.title('Пропущенные ошибки CRC (100 испытаний для каждого размера)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()