f1 = 17;
f2 = 17+4;
f3 = (17*2) + 1;

t = 0:0.001:1;

s1 = cos(2*pi*f1*t);
s2 = cos(2*pi*f2*t);
s3 = cos(2*pi*f3*t);

a = 4*s1 + 4*s2 + s3;
b = 3*s1 + s3;

% корреляция
corr1 = sum(s1 .* a);
corr2 = sum(s1 .* b);

% Нормализованная корреляция
norm1 = corr1 / sqrt(sum(s1.^2) * sum(a.^2));
norm2 = corr2 / sqrt(sum(s1.^2) * sum(b.^2));

fprintf('Корреляция s1-a: %.1f (норм: %.3f)\n', corr1, norm1);
fprintf('Корреляция s1-b: %.1f (норм: %.3f)\n', corr2, norm2);

a1 = [0.3 0.2 -0.1 4.2 -2 1.5 0];
b1 = [0.3 4 -2.2 1.6 0.1 0.1 0.2];

figure;
subplot(2,1,1);
plot(a1);
title('Массив a');
xlabel('Номер отсчета');
ylabel('Значение');
grid on;

subplot(2,1,2);
plot(b1);
title('Массив b');
xlabel('Номер отсчета');
ylabel('Значение');
grid on;

% взаимн корреляция
correlation = sum(a1 .* b1);
fprintf('Взаимная корреляция массивов a и b: %.2f\n', correlation);


results = [];
% Перебираем все сдвиги (от 0 до 6)
for shift = 0:6
    b_shift = circshift(b1, shift);  % circshift - циклически сдвигает массив b вправо
    corr = sum(a1 .* b_shift);% Считаем корреляцию: умножаем элементы a и b_shift и складываем
    results = [results, corr];
end

% Находим где была самая большая корреляция
[best_corr, best_shift] = max(results);
fprintf('Лучшая корреляция: %.2f при сдвиге %d\n', best_corr, best_shift-1);

figure;
plot(0:6, results, 'o-');
title('Корреляция от сдвига');
xlabel('Сдвиг'); ylabel('Корреляция');

b_best = circshift(b1, best_shift-1);
figure;
subplot(2,1,1); 
stem(a1, 'filled');
title('Массив a (максимальная корреляция)');
xlabel('Отсчет'); ylabel('Значение');
grid on;

subplot(2,1,2); 
stem(b_best, 'filled');
title(['Массив b (сдвиг ' num2str(best_shift-1) ', максимальная корреляция)']);
xlabel('Отсчет'); ylabel('Значение');
grid on;
