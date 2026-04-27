%% Первая последовательность Голда (x=12, y=19)
x_reg1 = [0 1 1 0 0]; % x = 12
y_reg1 = [1 0 0 1 1]; % y = 19

m_x1 = generate_m_seq(x_reg1, 'x');
m_y1 = generate_m_seq(y_reg1, 'y');
gold1 = xor(m_x1, m_y1);
L = length(gold1); % 31

fprintf('Первая последовательность Голда (x=12, y=19):\n');
fprintf('%d', gold1);
fprintf('\n\n');

fprintf('Сдвиг | бит 1 | бит 2 | бит 3 | Автокорреляция\n');
fprintf('------+-------+-------+-------+----------------\n');

R_auto_norm = zeros(1, L);
for tau = 0:L-1
    shifted = [gold1(L - tau + 1:end), gold1(1:L - tau)];
    matches = sum(gold1 == shifted);
    R_norm = (2 * matches - L) / L;
    R_auto_norm(tau + 1) = R_norm;
    fprintf('%5d | %5d | %5d | %5d | %+8.4f\n', ...
            tau, shifted(1), shifted(2), shifted(3), R_norm);
end

%% Вторая последовательность Голда (x=13, y=14)
x_reg2 = [0 1 1 0 1]; % x = 13
y_reg2 = [0 1 1 1 0]; % y = 14

m_x2 = generate_m_seq(x_reg2, 'x');
m_y2 = generate_m_seq(y_reg2, 'y');
gold2 = xor(m_x2, m_y2);

fprintf('\nВторая последовательность Голда (x=13, y=14):\n');
fprintf('%d', gold2);
fprintf('\n\n');

%% Нормированная взаимная корреляция
matches_cross = sum(gold1 == gold2);
R_cross_norm = (2 * matches_cross - L) / L;
fprintf('Взаимная корреляция (нормированная): %+8.4f\n', R_cross_norm);

fprintf('\nПроверка свойств\n');
check_balance_and_runs(m_x1, 'm-последовательность X1 (x=12)');
check_balance_and_runs(m_y1, 'm-последовательность Y1 (y=19)');
check_balance_and_runs(gold1, 'Последовательность Голда 1 (x=12, y=19)');

%% График
figure;
stem(0:L-1, R_auto_norm, 'filled');
title('Нормированная автокорреляционная функция последовательности Голда');
xlabel('Сдвиг \tau');
ylabel('R(\tau) (нормированная)');
grid on;
xlim([-1, L]);
ylim([-0.2, 1.1]);


%% Функция генерации m-последовательности
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

%%  проверка сбалансированности и циклов
function check_balance_and_runs(seq, name)
    L = length(seq);
    num_ones = sum(seq);
    num_zeros = L - num_ones;
    diff = abs(num_ones - num_zeros);
    
    % Анализ циклов (run-length) с учётом цикличности
    cyclic_seq = [seq, seq(1)]; 
    runs = [];
    current_len = 1;
    for i = 2:length(cyclic_seq)
        if cyclic_seq(i) == cyclic_seq(i-1)
            current_len = current_len + 1;
        else
            runs = [runs, current_len];
            current_len = 1;
        end
    end
    runs(runs > L) = []; 
    runs1 = sum(runs == 1);
    runs2 = sum(runs == 2);
    runs3 = sum(runs == 3);
    
    % Вывод
    fprintf('\n%s:\n', name);
    fprintf('  Баланс: %d единиц, %d нулей, разность = %d\n', num_ones, num_zeros, diff);
    if diff <= 1
        fprintf(' Сбалансированность выполнена\n');
    else
        fprintf(' Сбалансированность НЕ выполнена \n');
    end
    
    fprintf('  Циклы: длина 1 — %d, длина 2 — %d, длина 3 — %d\n', runs1, runs2, runs3);
    expected_runs1 = (L + 1) / 4; % для m-последовательности при L=31 → 8
    if abs(runs1 - expected_runs1) <= 1
        fprintf(' Распределение циклов корректно \n');
    else
        fprintf(' Распределение циклов отклоняется \n');
    end
end