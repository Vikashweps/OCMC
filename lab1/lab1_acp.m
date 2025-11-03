[y, Fs] = audioread('C:\Users\KL\Desktop\СибГУТИ\3 курс\ОСМС\lab1\voice.wav');

% Функция для квантования сигнала с визуализацией линий квантования
function [yq, l] = quantize_signal(y, bits, overload_factor)
    
    yo = y * overload_factor;% Усиливаем сигнал чтобы вызвать перегрузку
    y_clipped = max(-1, min(1, yo)); % Применяем ограничение (сигнал не может выйти за [-1, 1])
    
    % Определяем диапазон для квантования
    q_min = -1;
    q_max = 1;
    q_range = q_max - q_min;

    num_levels = 2^bits; % Количество уровней квантования
    delta = q_range / num_levels;% Шаг квантования
    l = q_min:delta:q_max;% Уровни квантования

    % Квантование сигнала
    y_normal = (y_clipped - q_min) / q_range;
    yq = round(y_normal * (num_levels - 1)) / (num_levels - 1);
    yq = yq * q_range + q_min;
    yq = yq * overload_factor;% Усиливаем обратно для отображения перегрузки
end

% Функция для вычисления коэффициентов Фурье 
function [an, bn] = fourier_coefficients(s, t, f, n_n)
    T = t(end) - t(1) + (t(2) - t(1)); % Полный период с учетом шага
    an = zeros(1, length(n_n));
    bn = zeros(1, length(n_n));
    for i = 1:length(n_n)
        n = n_n(i);
        sc = cos(2 * pi * n * f * t);  % опорное косинусное колебание
        ss = sin(2 * pi * n * f * t);
       
        s = s(:); % Убедимся, что все векторы - столбцы
        sc = sc(:);
        ss = ss(:);
        t = t(:);
        
        if n == 0
            an(i) = (1/T) * trapz(t, s .* sc);
            bn(i) = (1/T) * trapz(t, s .* ss);
        else
            an(i) = (2/T) * trapz(t, s .* sc);
            bn(i) = (2/T) * trapz(t, s .* ss);
        end
    end
    an = round(an, 4);
    bn = round(bn, 4);
    threshold = 1e-10;
    an(abs(an) < threshold) = 0.0;
    bn(abs(bn) < threshold) = 0.0;
end

% Функция для вычисления амплитудного и фазового спектра
function [An, fi] = compute_amplitude_phase(an, bn)
    An = sqrt(an.^2 + bn.^2);
    fi = atan2(-bn, an) * (180/pi);
    An = round(An, 4);
    fi = round(fi, 4);
end

% Функция для восстановления сигнала из коэффициентов Фурье
function s_reconstructed = reconstruct_signal(an, bn, t, f, n_n)
    s_reconstructed = zeros(size(t));
    
    for i = 1:length(n_n)
        n = n_n(i);
        if n == 0
            s_reconstructed = s_reconstructed + an(i);
        else
            s_reconstructed = s_reconstructed + an(i) * cos(2*pi*n*f*t) + bn(i) * sin(2*pi*n*f*t);
        end
    end
end

fprintf('Размер массива y: %d x %d\n', size(y, 1), size(y, 2));
y = y / max(abs(y)); % Нормализация к [-1, 1]

% Разрядности для анализа
bits_array = [3, 4, 5, 6]; % Разрядности АЦП
overload_factor =2.0; % Коэффициент перегрузки

% Найдем сегмент для анализа
segment_length = round(0.02 * Fs); % 20 мс
[max_val, max_idx] = max(abs(y));
start_idx = max(1, max_idx - round(segment_length/2));
end_idx = min(length(y), start_idx + segment_length - 1);

y_segment = y(start_idx:end_idx);
t_segment = (0:length(y_segment)-1) / Fs;

f0 = 100;% Основная частота и гармоники
n_n = 0:8;

% Анализ исходного сигнала
[an_orig, bn_orig] = fourier_coefficients(y_segment, t_segment, f0, n_n);
[An_orig, fi_orig] = compute_amplitude_phase(an_orig, bn_orig);

errors_original = zeros(size(bits_array));    % Ошибки для исходного сигнала
errors_overloaded = zeros(size(bits_array));  % Ошибки для перегруженного сигнала

for i = 1:length(bits_array)
    bits = bits_array(i);
    
    % Ошибка для исходного сигнала (без перегрузки)
    [y_quant_orig, ~] = quantize_signal(y_segment, bits, 1.0); % overload_factor = 1.0
    quant_er_orig = y_segment - y_quant_orig;
    er_orig(i) = mean(abs(quant_er_orig));
    
    % Ошибка для перегруженного сигнала
    y_q_over = quantize_signal(y_segment, bits, overload_factor);
    y_over = y_segment * overload_factor;
    y_q_over = y_over - y_q_over;
    er_over(i) = mean(abs(y_q_over));
end

figure('Name', 'Оцифровка исходного сигнала', 'Position', [100, 100, 1200, 800]);
% Исходный сигнал
subplot(2, 2, 1);
plot(t_segment, y_segment, 'k-', 'LineWidth', 2);
title('Исходный сигнал');
xlabel('Время, с'); ylabel('Амплитуда');
grid on;

% Оцифрованный сигнал для разных разрядностей с линиями квантования
colors = ['r', 'g', 'b', 'm'];
for i = 1:3 % Только 3,4,5 бит
    bits = bits_array(i);
    
    % Квантование БЕЗ перегрузки (только для визуализации)
    [y_quantized_no_overload, levels] = quantize_signal(y_segment, bits, 1.0); % overload_factor = 1.0
    
    subplot(2, 2, i+1);
    
    % Отображаем исходный сигнал
    plot(t_segment, y_segment, 'k-', 'LineWidth', 1); hold on;
    
    % Отображаем линии квантования
    for j = 1:length(levels)
        plot([t_segment(1), t_segment(end)], [levels(j), levels(j)], '--', ...
             'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.8);
    end
    
    % Отображаем квантованный сигнал
    plot(t_segment, y_quantized_no_overload, colors(i), 'LineWidth', 2);
    
    title(['Оцифровка: ' num2str(bits) ' бит (' num2str(2^bits) ' уровней)']);
    xlabel('Время, с'); ylabel('Амплитуда');
    legend('Исходный', 'Уровни квантования', 'Оцифрованный', 'Location', 'best');
    grid on;
end

figure('Name', 'Оцифровка перегруженного сигнала', 'Position', [100, 100, 1200, 800]);

% Исходный и перегруженный сигнал
subplot(2, 2, 1);
y_overloaded = y_segment * overload_factor;
plot(t_segment, y_segment, 'b-', 'LineWidth', 2); hold on;
plot(t_segment, y_overloaded, 'r-', 'LineWidth', 1.5);
title('Исходный и перегруженный сигнал');
xlabel('Время, с'); ylabel('Амплитуда');
legend('Исходный', ['Перегруженный ×' num2str(overload_factor)], 'Location', 'best');
grid on; 


% Оцифрованный перегруженный сигнал для разных разрядностей с линиями квантования
for i = 1:3 % Только 3,4,5 бит
    bits = bits_array(i);
    
    % Квантование С перегрузкой
    [y_quantized, levels] = quantize_signal(y_segment, bits, overload_factor);
    subplot(2, 2, i+1);
    % Отображаем перегруженный сигнал
    plot(t_segment, y_overloaded, 'k-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
    % Отображаем линии квантования 
    levels_scaled = levels * overload_factor;
    for j = 1:length(levels_scaled)
        plot([t_segment(1), t_segment(end)], [levels_scaled(j), levels_scaled(j)], '--', ...
             'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.8);
    end
    
    % Отображаем квантованный сигнал
    plot(t_segment, y_quantized, colors(i), 'LineWidth', 2);
    title(['Оцифровка перегруженного: ' num2str(bits) ' бит (' num2str(2^bits) ' уровней)']);
    xlabel('Время, с'); ylabel('Амплитуда');
    legend('Перегруженный', 'Уровни квантования', 'Оцифрованный', 'Location', 'best');
    grid on;
end

figure('Name', 'Влияние разрядности АЦП на спектр', 'Position', [100, 100, 1200, 800]);
for i = 1:length(bits_array)
    bits = bits_array(i);
    
    % Квантование исходного сигнала (без перегрузки)
    y_quantized_original = quantize_signal(y_segment, bits, 1.0); % overload_factor = 1.0
    
    % Расчет коэффициентов Фурье для квантованного исходного сигнала
    [an_quant_orig, bn_quant_orig] = fourier_coefficients(y_quantized_original, t_segment, f0, n_n);
    An_quant_orig = sqrt(an_quant_orig.^2 + bn_quant_orig.^2);
    
    subplot(2, 2, i);
    stem(n_n, An_orig, 'b-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    stem(n_n, An_quant_orig, 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
    title(['Разрядность: ' num2str(bits) ' бит (исходный сигнал)']);
    xlabel('n'); ylabel('An');
    legend('Исходный', 'Квантованный исходный', 'Location', 'best');
    grid on;
end