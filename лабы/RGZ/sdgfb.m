% Ввод данных
name = input('Имя латиницей: ', 's');
surname = input('Фамилия латиницей: ', 's');

text = [name ' ' surname];
fprintf('Исходный текст: "%s"\n', text);

% ASCII-кодер
bits = dec2bin(double(text), 8) - '0';

% Преобразование в одномерный массив
res = reshape(bits', 1, []);
L = length(res);

fprintf('Длина последовательности L = %d бит\n', L);
fprintf('Битовая последовательность:\n');
fprintf('%d', res);

% Параметры CRC
M = 7; % Длина CRC
G = [1 1 0 1 1 1 1 0]; % Порождающий полином

res_copy = [res zeros(1, M)]; % Добавляем M нулей

% Проходим по битам
for i = 1:length(res)
    if res_copy(i) == 1
        % XOR с полиномом
        for j = 0:7
            res_copy(i+j) = xor(res_copy(i+j), G(j+1));
        end
    end
end

% Берём последние M битов - это CRC
crc = res_copy(end-M+1:end);
fprintf('\nCRC-код (%d бит): ', M);
fprintf('%d', crc);

% Добавляем CRC к данным
res_crc = [res crc];

fprintf('\nДанные с CRC (%d бит):\n', length(res_crc));
fprintf('%d', res_crc);

% Проверка CRC 
check_res = res_crc; % Копируем для проверки
for i = 1:length(res_crc)-M  % ИСПРАВЛЕНО: правильный предел
    if check_res(i) == 1
        check_res(i:i+M) = xor(check_res(i:i+M), G);
    end
end

% Последние M битов должны быть нулями
crc_check = all(check_res(end-M+1:end) == 0);


% Поиск синхронизации с помощью корреляции
x_reg1 = [0 1 1 0 0]; % x = 12
y_reg1 = [1 0 0 1 1]; % y = 19

m_x1 = generate_m_seq(x_reg1, 'x');
m_y1 = generate_m_seq(y_reg1, 'y');
gold = xor(m_x1, m_y1); %последовательность Голда
G_len = length(gold); % 31
fprintf('\nПоследовательность Голда (%d бит):\n', G_len);
fprintf('%d', gold);

% АПСЕМПЛИНГ (должен быть ДО использования N!)
N = 5; % отсчетов на каждый бит
fprintf('\n\nN = %d отсчетов на бит\n', N);

gold_up = [];
for i = 1:length(gold)
    gold_up = [gold_up, gold(i)*ones(1, N)];
end
fprintf('Последовательность голда после апсемплинга(%d отсчетов):\n', length(gold_up));
fprintf('Первые 20: '); fprintf('%d', gold_up(1:20)); fprintf('\n');

tx_bits = [gold res_crc];

fprintf('\nПолный передаваемый сигнал (%d бит):\n', length(tx_bits));
fprintf('%d', tx_bits);

% АПСЕМПЛИНГ всего сигнала
up_signal = [];
for i = 1:length(tx_bits)
    up_signal = [up_signal, tx_bits(i)*ones(1, N)];
end
fprintf('\n\nАПСЕМПЛИНГ\n');
fprintf('Битов: %d\n', length(tx_bits));
fprintf('Отсчетов на бит: %d\n', N);
fprintf('Всего отсчетов: %d\n', length(up_signal));
fprintf('Первые 50 отсчетов:\n');
for i = 1:min(50, length(up_signal))
    fprintf('%d', up_signal(i));
    if mod(i, 5) == 0, fprintf(' '); end
    if mod(i, 25) == 0, fprintf('\n'); end
end

% ============= ПРОВЕРКА: БЕЗ ШУМА СНАЧАЛА =============
fprintf('\n\n=== ПРОВЕРКА БЕЗ ШУМА ===\n');
reserved_test = zeros(1, 2*N*(M+L+G_len));
position_test = 100; % Тестовая позиция
start_idx_test = position_test + 1;
end_idx_test = min(position_test + length(up_signal), length(reserved_test));

reserved_test(start_idx_test:end_idx_test) = up_signal(1:(end_idx_test-start_idx_test+1));

% Корреляция БЕЗ шума
Le_test = length(reserved_test);
corr_test = zeros(1, Le_test);

for shift = 0:Le_test-1
    start_idx_segment = shift + 1;
    end_idx_segment = shift + length(gold_up);
    
    if end_idx_segment <= Le_test
        segment = reserved_test(start_idx_segment:end_idx_segment);
        corr_test(shift+1) = sum(segment .* gold_up) / length(gold_up);
    else
        part1 = reserved_test(start_idx_segment:end);
        needed = length(gold_up) - length(part1);
        part2 = reserved_test(1:needed);
        segment = [part1, part2];
        corr_test(shift+1) = sum(segment .* gold_up) / length(gold_up);
    end
end

[corr_max_test, sync_sample_test] = max(corr_test);
fprintf('БЕЗ ШУМА: корреляция=%.4f, позиция=%d\n', corr_max_test, sync_sample_test);
fprintf('Ожидаемая позиция: %d\n', start_idx_test);

if abs(corr_max_test - 1.0) < 0.01
    fprintf('✓ Корреляция идеальная без шума\n');
else
    fprintf('⚠ Проблема: даже без шума корреляция %.4f (должна быть 1.0)\n', corr_max_test);
end
