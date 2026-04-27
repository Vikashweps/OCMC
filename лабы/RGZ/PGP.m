name = input('Имя латиницей: ', 's');
surname = input('Фамилия латиницей: ', 's');
text = [name ' ' surname];
fprintf('Исходный текст: "%s"\n', text);

%% ASCII-кодер
bits = dec2bin(double(text), 8) - '0';
res = reshape(bits', 1, []); % Преобразование в одномерный массив
L = length(res); % L-длина последовательности
fprintf('Длина последовательности L = %d бит\n', L);
fprintf('%d', res);

%% CRC
M = 7; % Длина CRC
poli = [1 1 0 1 1 1 1 0]; % Порождающий полином
res_copy = [res zeros(1, M)]; % Добавляем M нулей

% Проходим по битам
for i = 1:length(res)
    if res_copy(i) == 1
        % XOR с полиномом
        for j = 0:7
            res_copy(i+j) = xor(res_copy(i+j), poli(j+1));
        end
    end
end
crc = res_copy(end-M+1:end);% Берём последние M битов - это CRC
fprintf('\nCRC-код (%d бит): ', M);
fprintf('%d', crc);
res_crc = [res crc]; % Добавляем CRC к данным

fprintf('\nДанные с CRC (%d бит):\n', length(res_crc));
fprintf('%d', res_crc);

%% последовательность Голда
x_reg = [0 1 1 0 0]; % x = 12
y_reg = [1 0 0 1 1]; % y = 19
m_x = generate_m_seq(x_reg, 'x');
m_y = generate_m_seq(y_reg, 'y');
gold = xor(m_x, m_y);
G = length(gold); % 31

fprintf('\nПоследовательность Голда (%d бит):\n', G);
fprintf('%d', gold);

tx_bits = [gold res_crc];
fprintf('\nПолный передаваемый сигнал (%d бит):\n', length(tx_bits));
fprintf('%d', tx_bits);

%% апсемплинг
N = 10; % отсчетов на бит
up_signal = repelem(tx_bits, N);
gold_up = repelem(gold, N);

fprintf('\n\nАПСЕМПЛИНГ\n', N);

fprintf('Всего отсчетов: %d\n', length(up_signal));
fprintf('Первые 50 отсчетов:\n');
for i = 1:min(50, length(up_signal))
    fprintf('%d', up_signal(i));
    if mod(i, 5) == 0, fprintf(' '); end
    if mod(i, 25) == 0, fprintf('\n'); end
end

fprintf('\nПоследовательность голда после апсемплинга(%d бит):\n', length(gold_up));
fprintf('%d', gold_up);

%% создаем сторой массив и вставка
reserved = zeros(1, 2*N*(M+L+G));
fprintf('\nДлина второго массива: %d бит',length(reserved));
position = input(sprintf('\nВведите позицию для вставки (0-%d): ', (N*(M+L+G))));
start_idx = position + 1;

% Проверяем, чтобы не выйти за границы массива
if start_idx <= length(up_signal)
    end_idx = min(position + length(up_signal), length(reserved));
    reserved(start_idx:end_idx) = up_signal(1:(end_idx-start_idx+1));
else
    fprintf('Ошибка: позиция выходит за границы массива\n');
end

%% Шум
sigma = input('σ (float): ');
noise = sigma * randn(1, length(reserved)); % Генерируем шум по нормальному распределению
fprintf('Размер массива шума: %d отсчетов\n', length(noise));
received_signal = reserved + noise; %  Добавляем шум к сигналу

%% корреляция
Le = length(received_signal);
Lg = length(gold_up);
norm_corr = zeros(1, Le);

for shift = 0:Le-1
    % Берем сегмент сигнала длиной Lg, начиная с позиции shift+1
    start_idx_segment = shift + 1;
    end_idx_segment = shift + Lg;
    
    if end_idx_segment <= Le
        % Если сегмент полностью внутри сигнала
        segment = received_signal(start_idx_segment:end_idx_segment);
        norm_corr(shift+1) = sum(segment .* gold_up) / sqrt(sum(segment.^2) * sum(gold_up.^2));
    else
        % Если сегмент выходит за границы - берем циклически
        part1 = received_signal(start_idx_segment:end);
        needed = Lg - length(part1);
        part2 = received_signal(1:needed);
        segment = [part1, part2];
        norm_corr(shift+1) = sum(segment .* gold_up) / sqrt(sum(segment.^2) * sum(gold_up.^2));
    end
end

[corr_max, sync_sample] = max(norm_corr);
fprintf('корреляция=%.4f, позиция=%d\n', corr_max, sync_sample);

%% Удаляем отсчеты до синхры
signal_clean = received_signal(sync_sample:end);
fprintf('Удалено %d отсчетов\n', sync_sample-1);

%% ДЕМОДУЛЯЦИЯ
P = 0.5;
total_bits = G + L + M;
fprintf('Всего битов для демодуляции: %d\n', total_bits);

bits_out = zeros(1, total_bits);

% Проверяем хватает ли отсчетов
if length(signal_clean) < total_bits * N
    fprintf('ОШИБКА: Недостаточно отсчетов!\n');
    fprintf('Есть: %d, нужно: %d\n', length(signal_clean), total_bits*N);
else
    for i = 1:total_bits
        start = (i-1)*10 + 1;
        stop = start + 9;
        five_samples = signal_clean(start:stop);
        avg = mean(five_samples);
        bits_out(i) = avg > P;
    end
end
    fprintf('Принято битов (%d):\n', length(bits_out));
    fprintf('%d', bits_out);
    fprintf('\n');
    
    %% Разделяем на части
    data_crc_rx = bits_out(G+1:end); % данные без последовательности синхронизации
    crc_rx = data_crc_rx(end-6:end); %  CRC
    data_rx = data_crc_rx(1:end-7);  % данные (без CRC)
    fprintf('Данные+CRC: ');
    fprintf('%d', data_crc_rx);
    fprintf('\n');
    
  %% Проверка CRC
  check_data = [data_rx zeros(1,7)];

% Вычисляем CRC для принятых данных
for i = 1:length(data_rx)
    if check_data(i) == 1
        check_data(i:i+7) = xor(check_data(i:i+7), poli);
    end
end

% Получаем вычисленный CRC
calculated_crc = check_data(end-6:end);

% Сравниваем с принятым CRC
if all(calculated_crc == crc_rx)
    fprintf('CRC верный\n');
    fprintf('Принятый CRC: \n'); fprintf(' %d', crc_rx);
    fprintf('\nВычисленный CRC: \n' ); fprintf(' %d', calculated_crc);
else
    fprintf('Ошибка CRC\n');
    fprintf('Принятый CRC: \n'); fprintf(' %d', crc_rx);
    fprintf('\nВычисленный CRC: \n' ); fprintf(' %d', calculated_crc);
end

% Декодирование для проверки
fprintf('\nПроверка декодированием: "%s"\n', char(bin2dec(reshape(char(data_rx' + '0'), 8, [])')));

% Графики
figure();
stairs(res, 'LineWidth', 2);
title('Представление битовой последовательности');
xlabel('Номер бита');
ylabel('Значение бита');
ylim([-0.1 1.1]);
grid on;

figure();
stairs(tx_bits, 'b', 'LineWidth', 1.5);
hold on;
plot([31.5 31.5], [-0.1 1.1], 'r--', 'LineWidth', 2);
title('Сигнал: Голд + данные + CRC');
xlabel('Бит'); ylabel('Значение');
ylim([-0.1 1.1]); grid on;
legend('Биты', 'Граница синхр.', 'Location', 'best');

% Временной сигнал 
figure();
stairs(up_signal, 'LineWidth', 1);
title(sprintf('Временной сигнал после апсемплинга (N=%d, всего %d отсчетов)', N, length(up_signal)));
xlabel('Номер отсчета');
ylabel('Амплитуда');
ylim([-0.1 1.1]);
xlim([0, (length(up_signal) + 15)]);
grid on;

figure();
stairs(reserved, 'b-', 'LineWidth', 1);
hold on;
if start_idx <= end_idx 
    plot(start_idx:end_idx, reserved(start_idx:end_idx), 'r-', 'LineWidth', 2);
end
title('Массив с вставленным сигналом');
xlabel('Отсчет'); ylabel('Амплитуда');
ylim([-0.1 1.1]);
xlim([0, (length(reserved) + 15)]);
grid on;

figure();
plot(received_signal, 'g-', 'LineWidth', 1);
hold on;
title('Принятый зашумленный сигнал');
xlabel('Номер отсчета');
ylabel('Амплитуда');
ylim([min(received_signal)-0.2, max(received_signal)+0.2]);
grid on;

% график корреляционной функции
figure();
plot(norm_corr, 'b-', 'LineWidth', 1);
hold on;
plot([sync_sample, sync_sample], [min(norm_corr), max(norm_corr)], 'r--', 'LineWidth', 2);
title('Нормированная корреляционная функция');
xlabel('Позиция в сигнале');
ylabel('Корреляция');
grid on;

%% СПЕКТРЫ 
figure;

% Спектр передаваемого
subplot(2,1,1);
spectrum_tx = abs(fft(up_signal));
n_tx = 1:length(spectrum_tx)/2;
stem(n_tx, spectrum_tx(1:end/2));
title('Спектр передаваемого сигнала');
xlabel('Номер гармоники');
ylabel('Амплитуда');
xlim([-0.1 101]);
grid on;

% Спектр принимаемого
subplot(2,1,2);
spectrum_rx = abs(fft(received_signal));
n_rx = 1:length(spectrum_rx)/2;
stem(n_rx, spectrum_rx(1:end/2));
xlim([-0.1 101]);
title('Спектр принимаемого сигнала');
xlabel('Номер гармоники');
ylabel('Амплитуда');
grid on;

%% Генерация m-последовательности
function seq = generate_m_seq(state, type)
    n = length(state);
    seq = zeros(1, 2^n - 1);
    reg = state(:)';

    for i = 1:(2^n - 1)
        seq(i) = reg(end);
        if strcmp(type, 'x')
            new_bit = xor(reg(3), reg(5)); % x4 XOR x5
        else
            new_bit = xor(reg(3), reg(5)); % y3 XOR y5
        end
        reg = [new_bit, reg(1:end-1)];
    end
end