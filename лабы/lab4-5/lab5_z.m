data_sizes = 1000:1000:10000;   % 5 точек: 1k, 3k, ..., 10k (для скорости)
crc_lengths = 3:8;
num_trials = 200;

% ПЛОХИЕ полиномы: все заканчиваются на 0 (нет x^0 = 1)
polynomials = {
    [1 1 1 0], ...        % CRC-3: x^3 + x^2 + x          (без +1)
    [1 1 0 1 0], ...      % CRC-4: x^4 + x^3 + x          (без +1)
    [1 1 0 1 0 0], ...    % CRC-5: x^5 + x^4 + x^2        (без +1)
    [1 1 0 0 0 1 0], ...  % CRC-6: x^6 + x^5 + x          (без +1)
    [1 1 0 0 1 0 0 0], ...% CRC-7: x^7 + x^6 + x^3        (без +1)
    [1 1 0 0 0 0 1 1 0]   % CRC-8: x^8 + x^7 + x^2 + x    (без +1)
};

detection_rates = zeros(length(crc_lengths), length(data_sizes));

for ci = 1:length(crc_lengths)
    n = crc_lengths(ci);
    poly = polynomials{ci};
    fprintf('Тестируем CRC-%d (плохой полином)\n', n);
    
    for di = 1:length(data_sizes)
        L = data_sizes(di);
        total_bits = L + n;
        detected = 0;
        
        for t = 1:num_trials
            % Генерация данных
            data = randi([0 1], 1, L);
            
            % Вычисление CRC
            crc = calculateCRC(data, poly);
            
            % Формирование пакета
            packet = [data, crc];
            
            % ОДИНОЧНАЯ ОШИБКА
            pos = randi([1, total_bits]);
            packet_err = packet;
            packet_err(pos) = 1 - packet_err(pos);
            
            % Проверка
            if ~checkCRC(packet_err, poly)
                detected = detected + 1;
            end
        end
        
        detection_rates(ci, di) = detected / num_trials;
    end
end

%% Построение графика
figure('Position', [100, 100, 800, 500]);
colors = lines(length(crc_lengths));

for ci = 1:length(crc_lengths)
    plot(data_sizes, detection_rates(ci, :), 'Color', colors(ci,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('CRC-%d', crc_lengths(ci)));
    hold on;
end

xlabel('Размер полезной нагрузки (бит)');
ylabel('Вероятность обнаружения ОДИНОЧНОЙ ошибки');
title('Эффективность CRC с ПЛОХИМИ полиномами (без x^0)');
grid on;
legend('Location', 'southeast');
ylim([0.7, 1.0]);  % фокус на области, где есть различия

function crc = calculateCRC(data, poly)
    n = length(poly) - 1;
    crc_reg = zeros(1, n);
    for i = 1:length(data)
        fb = xor(crc_reg(1), data(i));
        crc_reg(1:end-1) = crc_reg(2:end);
        crc_reg(end) = 0;
        if fb
            crc_reg = xor(crc_reg, poly(2:end));
        end
    end
    crc = crc_reg;
end

function ok = checkCRC(packet, poly)
    n = length(poly) - 1;
    crc_reg = zeros(1, n);
    for i = 1:length(packet)
        fb = xor(crc_reg(1), packet(i));
        crc_reg(1:end-1) = crc_reg(2:end);
        crc_reg(end) = 0;
        if fb
            crc_reg = xor(crc_reg, poly(2:end));
        end
    end
    ok = all(crc_reg == 0);
end