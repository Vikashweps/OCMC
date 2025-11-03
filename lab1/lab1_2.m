[y, Fs] = audioread('C:\Users\KL\Desktop\СибГУТИ\3 курс\ОСМС\lab1\voice.wav');

if exist('voice.wav', 'file') == 2
    fprintf('Файл найден!\n');
    [y, Fs] = audioread('voice.wav');
else
    fprintf('Файл не найден в текущей папке\n');
    fprintf('Текущая папка: %s\n', pwd);
end
fprintf('Размер массива y: %d x %d\n', size(y, 1), size(y, 2));
fprintf('Частота дискретизации Fs: %d Гц\n', Fs);

N = length(y);
duration = N / Fs;
Fs_calc = N / duration;

fprintf('Число отсчётов (N): %d\n', N);
fprintf('Длительность записи: %.3f сек\n', duration);
fprintf('Частота дискретизации (расчётная): %.1f Гц\n', Fs_calc);

y1 = downsample(y, 10);
Fs_new = Fs / 10;

fprintf('\nИсходная частота дискретизации: %d Гц\n', Fs);
fprintf('Коэффициент прореживания: 10\n');
fprintf('Новая частота дискретизации: %d Гц\n', Fs_new);
fprintf('Исходное количество отсчетов: %d\n', length(y));
fprintf('После прореживания: %d отсчетов\n', length(y1));

fprintf('\nВоспроизведение прореженного сигнала...\n');
audio = audioplayer(y1, Fs_new);
play(audio);

Y_orig = fft(y);
L_orig = length(Y_orig);
P2_orig = abs(Y_orig/L_orig);
P1_orig = P2_orig(1:L_orig/2+1);
P1_orig(2:end-1) = 2*P1_orig(2:end-1);
f_orig = Fs*(0:(L_orig/2))/L_orig;

Y_down = fft(y1);
L_down = length(Y_down);
P2_down = abs(Y_down/L_down);
P1_down = P2_down(1:L_down/2+1);
P1_down(2:end-1) = 2*P1_down(2:end-1);
f_down = Fs_new*(0:(L_down/2))/L_down;

P1_orig_db = 20*log10(P1_orig + eps);
P1_down_db = 20*log10(P1_down + eps);

threshold_orig_db = max(P1_orig_db) - 40;
mask_orig = P1_orig_db > threshold_orig_db;
f_min_orig = f_orig(find(mask_orig, 1, 'first'));
f_max_orig = f_orig(find(mask_orig, 1, 'last'));
bandwidth_orig = f_max_orig - f_min_orig;

threshold_down_db = max(P1_down_db) - 40;
mask_down = P1_down_db > threshold_down_db;
f_min_down = f_down(find(mask_down, 1, 'first'));
f_max_down = f_down(find(mask_down, 1, 'last'));
bandwidth_down = f_max_down - f_min_down;

figure('Name', '12. Сравнение амплитудных спектров', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(f_orig, P1_orig, 'b', 'LineWidth', 1);
title(['Оригинальный сигнал (Fs = ' num2str(Fs) ' Гц)']);
xlabel('Частота, Гц'); ylabel('Амплитуда');
xlim([0, 8000]);
grid on;
hold on;
plot([f_min_orig, f_max_orig], [max(P1_orig)*0.5, max(P1_orig)*0.5], 'r--', 'LineWidth', 2);
legend('Спектр', 'Ширина спектра');

subplot(2,2,2);
plot(f_down, P1_down, 'r', 'LineWidth', 1);
title(['Прореженный сигнал (Fs = ' num2str(Fs_new) ' Гц)']);
xlabel('Частота, Гц'); ylabel('Амплитуда');
xlim([0, Fs_new/2]);
grid on;
hold on;
plot([f_min_down, f_max_down], [max(P1_down)*0.5, max(P1_down)*0.5], 'b--', 'LineWidth', 2);
legend('Спектр', 'Ширина спектра');

subplot(2,2,[3,4]);
plot(f_orig, P1_orig, 'b', 'LineWidth', 1.5); hold on;
plot(f_down, P1_down, 'r', 'LineWidth', 1);
title('Сравнение спектров оригинального и прореженного сигналов');
xlabel('Частота, Гц'); ylabel('Амплитуда');
xlim([0, 8000]);
legend(['Оригинал (Fs=' num2str(Fs) 'Гц)'], ['Прореженный (Fs=' num2str(Fs_new) 'Гц)']);
grid on;

freq_ticks = [50, 60, 80, 100, 120, 150, 200, 240, 300, 400, 500, 700, 1000, 1300, 1600, 2000, 3000, 4000, 6000, 8000, 9000];
db_ticks = -90:5:-25;

figure('Name', 'график частотной характеристики', 'Position', [100, 100, 1400, 700]);

semilogx(f_orig, P1_orig_db, 'b-', 'LineWidth', 2);
hold on;
xlim([50, 9000]);
ylim([-90, -25]);
xticks(freq_ticks);
yticks(db_ticks);

ax = gca;
ax.XTickLabelRotation = 45;
grid on;
grid minor;
xlabel('Частота, Гц', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Уровень, дБ', 'FontSize', 12, 'FontWeight', 'bold');
title('Плотный график частотной характеристики', 'FontSize', 14, 'FontWeight', 'bold');
legend(['Оригинал (Fs=' num2str(Fs) 'Гц)'], 'Location', 'southwest', 'FontSize', 10);
hold off;

figure('Name', 'Сравнение спектров в логарифмическом масштабе', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
semilogx(f_orig, P1_orig_db, 'b', 'LineWidth', 1);
title(['Оригинальный сигнал (Fs = ' num2str(Fs) ' Гц) в дБ']);
xlabel('Частота, Гц'); ylabel('Амплитуда, дБ');
xlim([50, 8000]);
ylim([-85, -30]);
xticks(freq_ticks(freq_ticks <= 8000));
grid on;

subplot(2,2,2);
semilogx(f_down, P1_down_db, 'r', 'LineWidth', 1);
title(['Прореженный сигнал (Fs = ' num2str(Fs_new) ' Гц) в дБ']);
xlabel('Частота, Гц'); ylabel('Амплитуда, дБ');
xlim([50, Fs_new/2]);
ylim([-85, -30]);
xticks(freq_ticks(freq_ticks <= Fs_new/2));
grid on;

subplot(2,2,[3,4]);
semilogx(f_orig, P1_orig_db, 'b', 'LineWidth', 1.5); hold on;
semilogx(f_down, P1_down_db, 'r', 'LineWidth', 1);
title('Сравнение спектров в логарифмическом масштабе (дБ)');
xlabel('Частота, Гц'); ylabel('Амплитуда, дБ');
xlim([50, 8000]);
ylim([-85, -30]);
xticks(freq_ticks(freq_ticks <= 8000));
legend(['Оригинал (Fs=' num2str(Fs) 'Гц)'], ['Прореженный (Fs=' num2str(Fs_new) 'Гц)'], 'Location', 'best');
grid on;