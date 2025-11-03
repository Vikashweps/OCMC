f = 13;  % Гц

%% непрерывный сигнал
t = 0:0.0001:1;  
y = sin(2*pi*f*t) + sin(10*pi*f*t);

figure('Name', '1. Непрерывный сигнал');
plot(t, y, 'LineWidth', 1.5);
title(['Сигнал y(t) = sin(2\pi \cdot ', num2str(f), ' \cdot t) + sin(10\pi \cdot ', num2str(f), ' \cdot t)']);
xlabel('Время, с');
ylabel('Амплитуда');
grid on;
xlim([0, 0.2]);  

%% Макс частота и частота дискретизации по Котельникову
f1 = f;          % 13 Гц
f2 = 5 * f;      % 65 Гц (т.к. sin(10πft) = sin(2π·5f·t))
f_max = max(f1, f2);  % 65 Гц
fd_min = 2 * f_max;   % 130 Гц

fprintf('Максимальная частота в спектре: %d Гц\n', f_max);
fprintf('Минимальная частота дискретизации (Котельников): %d Гц\n', fd_min);

%% Оцифровка сигнала (дискретные отсчёты) 
fd = fd_min +1;
T = 1;           % 1 секунда
N = fd * T;  % 131 отсчёт

td = (0:N-1) / fd; 
yd = sin(2*pi*f*td) + sin(10*pi*f*td);  

figure('Name', '4. Дискретные отсчёты (fs = 131 Гц)');
stem(td, yd, 'filled');
title(['Дискретные отсчёты (f_s = ', num2str(fs), ' Гц)']);
xlabel('Время, с');
ylabel('Амплитуда');
grid on;
xlim([0, 0.2]);

%%  Прямое ДПФ (FFT), ширина спектра и объём памяти
Y = fft(yd);

f_axis = (0:N/2) * fd / N;
Y_mag = abs(Y(1:N/2+1)) / N * 2;
Y_mag(1) = Y_mag(1) / 2;  % постоянная составляющая

spectrum_width = f2 - f1;  % 52 Гц

memory_64 = N * 8;   % байт

fprintf('\n Результаты ПДПФ при f_s = %d Гц \n', fs);
fprintf('Ширина спектра: %d Гц (от %d до %d Гц)\n', spectrum_width, f_min, f_max);
fprintf('Объём памяти (float64): %d байт\n', memory_64);

figure('Name', '5. Амплитудный спектр (fs = 131 Гц)');
stem(f_axis, Y_mag, 'filled');
title(['Амплитудный спектр сигнала (f_s = ', num2str(fs), ' Гц)']);
xlabel('Частота, Гц');
ylabel('Амплитуда');
grid on;
xlim([0, 80]);
ylim([0, 1.2]);

%% Восстановление сигнала 
t_interp = 0:0.0001:0.2;
y_interp = interp1(td, yd, t_interp, 'linear');

y_orig = sin(2*pi*f*t_interp) + sin(10*pi*f*t_interp);

figure('Name', '6. Восстановление (fs = 131 Гц)');
plot(t_interp, y_orig, 'b-', 'LineWidth', 1.5); hold on;
plot(t_interp, y_interp, 'r--', 'LineWidth', 1.2);
stem(td(td <= 0.2), yd(td <= 0.2), 'k', 'filled');
legend('Оригинальный сигнал', 'Восстановленный (линейная)', 'Отсчёты');
xlabel('Время, с');
ylabel('Амплитуда');
title(['Сравнение оригинала и восстановленного (f_s = ', num2str(fs), ' Гц)']);
grid on;
xlim([0, 0.2]);
%% 7пункт
fd2 = fd_min*4;
N2 = round(fd2 *T);
t2 = (0:N2-1) / fd2;
y2 = sin(2*pi*f*t2) + sin(10*pi*f*t2);

figure('Name', '4. Дискретные отсчёты (fs = 520 Гц)');
stem(t2, y2, 'filled');
title(['Дискретные отсчёты (f_s = ', num2str(fd2), ' Гц)']);
xlabel('Время, с');
ylabel('Амплитуда');
grid on;
xlim([0, 0.2]);

% ДПФ для fd2
Y2 = fft(y2);
f_axis2 = (0:N2/2) * fd2 / N2;
Y_mag2 = abs(Y2(1:N2/2+1)) / N2 * 2;
Y_mag2(1) = Y_mag2(1) / 2;

% Восстановление
y2_interp = interp1(t2, y2, t_interp, 'linear');

% Вывод
memory_64_2 = N2 * 8;
fprintf('\n Результаты при f_d = %d Гц (в 4 раза выше) \n', fd2);
fprintf('Ширина спектра: %d Гц (от %d до %d Гц)\n', spectrum_width, f_min, f_max);
fprintf('Объём памяти (float64): %d байт\n', memory_64_2);

% График спектра для fd2
figure('Name', '7a. Спектр при fd = 520 Гц');
stem(f_axis2, Y_mag2, 'filled');
title(['Амплитудный спектр (f_d = ', num2str(fd2), ' Гц)']);
xlabel('Частота, Гц');
ylabel('Амплитуда');
grid on;
xlim([0, 80]);
ylim([0, 1.2]);

% График восстановления для fd2
figure('Name', '7b. Восстановление при fd = 520 Гц');
plot(t_interp, y_orig, 'b-', 'LineWidth', 1.5); hold on;
plot(t_interp, y2_interp, 'r--', 'LineWidth', 1.2);
stem(t2(t2 <= 0.2), y2(t2 <= 0.2), 'k', 'filled');
legend('Оригинальный сигнал', 'Восстановленный (линейная)', 'Отсчёты');
xlabel('Время, с');
ylabel('Амплитуда');
title(['Восстановление (f_d = ', num2str(fd2), ' Гц)']);
grid on;
xlim([0, 0.2]);
