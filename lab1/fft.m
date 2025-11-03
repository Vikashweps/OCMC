f = 5;
Ts = 0.001;
N_values = [128, 256, 512, 1024, 2048, 4096];  % степени двойки
num_iter = 500;  % количество итераций для усреднения

time_my = zeros(size(N_values));
time_fft = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    t = (-N/2 * Ts) : Ts : (N/2 * Ts - Ts);
    s = 2 * cos(2 * pi * f * t);

    total_my = 0;
    total_fft = 0;

    for k = 1:num_iter
        %  Прямой метод 
        tic;
        T = 1 / f;
        n_n = 0:4;
        for n = n_n
            sc = cos(2 * pi * n * f * t);
            ss = sin(2 * pi * n * f * t);
            if n == 0
                a0 = (1 / T) * trapz(t, s .* sc);
                b0 = (1 / T) * trapz(t, s .* ss);
            else
                a = (2 / T) * trapz(t, s .* sc);
                b = (2 / T) * trapz(t, s .* ss);
            end
        end
        total_my = total_my + toc;

        % FFT 
        tic;
        fft(s);
        total_fft = total_fft + toc;
    end

    time_my(i) = total_my / num_iter;
    time_fft(i) = total_fft / num_iter;
end

figure;
semilogy(N_values, time_my, 'o-', 'DisplayName', 'Прямой метод');
hold on;
semilogy(N_values, time_fft, 's-', 'DisplayName', 'БПФ (FFT)');
xlabel('Число точек N');
ylabel('Среднее время выполнения (сек), лог. шкала');
title(['Сравнение производительности (Ts = ', num2str(Ts), ' с, ', num2str(num_iter), ' итераций)']);
legend('Location', 'best');
grid on;
set(gca, 'XTick', N_values);