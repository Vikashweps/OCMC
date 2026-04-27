a = [3 2 -1 4 2 -2 1 0 5 1 3 2];
results1 = [];
results2 = [];
lags = 0:length(a)-1;

for shift = 0:length(a)-1
    % Циклический сдвиг
    a_shift = zeros(1, length(a));
    for i = 1:length(a)
        new_pos = mod(i - 1 - shift, length(a)) + 1;
        a_shift(i) = a(new_pos);
    end
    
    % xcorr для циклически сдвинутого массива
    [c_temp, lags_temp] = xcorr(a, a_shift, 'normalized');
    zero_lag_index = find(lags_temp == 0, 1);
    results1 = [results1, c_temp(zero_lag_index)];
    
    % ручной расчет
    corr = sum(a .* a_shift);
    norm_corr = corr / (sqrt(sum(a.^2)) * sqrt(sum(a_shift.^2)));
    results2 = [results2, norm_corr];
end


fprintf('Сдвиг | xcorr метод | Ручной метод\n');
fprintf('------|-------------|-------------\n');
for i = 1:length(a)
    fprintf('  %2d  |    %8.4f   |   %8.4f\n', i-1, results1(i), results2(i));
end

% Графики
figure;
subplot(2,1,1);
stem(lags, results1, 'filled');
title('Метод с xcorr');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on;

subplot(2,1,2);
stem(lags, results2, 'filled', 'r');
title('Ручной метод');
xlabel('Сдвиг'); ylabel('Корреляция');
grid on;