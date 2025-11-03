TxPowerBS = 46;          % дБм, мощность передатчика BS
TxPowerUE = 24;          % дБм, мощность передатчика UE
AntGainBS = 21;          % дБи, коэффициент усиления антенны BS
PenetrationM = 15;       % дБ, запас на проникновение
IM = 1;                  % дБ, запас на интерференцию
f = 1.8;                % ГГц, частота
BW_UL = 10e6;            % Гц, полоса UL
BW_DL = 20e6;            % Гц, полоса DL
NF_BS = 2.4;             % дБ, коэффициент шума BS
NF_UE = 6;               % дБ, коэффициент шума UE
SINR_DL = 2;             % дБ, требуемое SINR для DL
SINR_UL = 4;             % дБ, требуемое SINR для UL
MIMO = 2;       % количество антенн MIMO
MIMOGain = 3;    % MIMO = 2 ант усилит сигнал на 3дБ или в 2 раза
S = 100;         % кв.км, общая площадь
S_trk = 4;          % кв.км, площадь помещений

% Тепловой шум
TN_UL = -174 + 10*log10(BW_UL);
TN_DL = -174 + 10*log10(BW_DL);
% (UL)
RxSensBS = NF_BS + TN_UL + SINR_UL;
% (DL)  
RxSensUE = NF_UE + TN_DL + SINR_DL;

feederLoss = 0.5 + 0.4 + 2;

% Бюджет DL
MAPL_DL = TxPowerBS - feederLoss + AntGainBS + MIMOGain - IM - PenetrationM - RxSensUE;

% Бюджет UL
MAPL_UL = TxPowerUE - feederLoss + AntGainBS + MIMOGain - IM - PenetrationM - RxSensBS;

fprintf('MAPL_DL = %.2f дБ\nMAPL_UL = %.2f дБ\n', MAPL_DL, MAPL_UL);

% Диапазон расстояний
d_m = [200, 500, 1000, 2000, 5000]; % расстояния в метрах
d_km = d_m / 1000;

% Параметры для моделей
h_BS = 30;    % м, высота антенны BS
h_MS = 1.5;   % м, высота антенны MS

% Коэффициенты для диапазона 1500-2000 МГц
A = 46.3;
B = 33.9;

% UMiNLOS (микросоты)
PL_umi = 26*log10(f) + 22.7 + 36.7*log10(d_m); % f в ГГц, d в м

%COST231
a_hMS = 3.2*(log10(11.75*h_MS))^2 - 4.97;% для городской местности
L_clutter = 0;      % для городской местности
s = 44.9 - 6.55*log10(h_BS);    % для d >= 1км

PL_cost231 = A + B*log10(f*1000) - 13.82*log10(h_BS) - a_hMS + s*log10(d_km) + L_clutter;

% Walfish-Ikegami
PL_WI = 42.6 + 20 * log10(f*1000) + 26 * log10(d_km);

r_UL = 550;
r_DL = 1740;
R = min(r_UL, r_DL); % в м
S_site = 1.95 * (R/1000)^2;
site_trk = ceil(S_trk / S_site);
site_o = ceil(S / S_site);

fprintf('Радиус соты: %.1f м\n', R);
fprintf('Площадь сайта: %.3f км²\n', S_site);
fprintf('Общее число БС: %d\n', site_o);
fprintf('Число БС для микрозоны: %d\n', site_trk);
%% Вывод результатов в таблицу
fprintf('Зависимость потерь радиосигнала от расстояния\n');
fprintf('Расстояние | UMiNLOS | COST231 Hata | Walfish-Ikegami\n');
fprintf('   (м)     |  (дБ)   |    (дБ)      |      (дБ)     \n');
fprintf('-----------|---------|--------------|----------------\n');

for i = 1:length(d_m)
    fprintf('   %4d    | %7.1f | %12.1f | %14.1f\n', ...
            d_m(i), PL_umi(i), PL_cost231(i), PL_WI(i));
end

% Построение графика
figure('Position', [200, 200, 1000, 600]);

% График потерь с радиусами

plot(d_m, PL_umi, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'UMiNLOS');
hold on;
plot(d_m, PL_cost231, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'COST231 Hata');
plot(d_m, PL_WI, 'g-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Walfish-Ikegami');

% Линии MAPL
yline(MAPL_UL, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('MAPL_{UL} = %.1f дБ', MAPL_UL));
yline(MAPL_DL, 'k:', 'LineWidth', 1.5, 'DisplayName', sprintf('MAPL_{DL} = %.1f дБ', MAPL_DL));
grid on;
xlabel('Расстояние, м');
ylabel('Потери сигнала, дБ');
title('Зависимость потерь от расстояния');
legend('Location', 'northwest');


% Параметры для графика
frequencies = 400:100:3000;  % МГц
heights = 0:10:200;          % высота BS в метрах
h_MS = 1.5;                  % высота пользователя

% Массивы для радиусов
R_DL = zeros(length(frequencies), length(heights));
R_UL = zeros(length(frequencies), length(heights));

% Расчет радиусов
for i = 1:length(frequencies)
    for j = 1:length(heights)
        f_ghz = frequencies(i) / 1000;  % в ГГц
        h_BS = heights(j);
        
        % Параметры модели затухания
        a_hMS = 3.2*(log10(11.75*h_MS))^2 - 4.97;
        s = 44.9 - 6.55*log10(h_BS);
        A = 46.3;
        B = 33.9;
        
        % Расчет радиуса DL
        PL_const_DL = A + B*log10(f_ghz*1000) - 13.82*log10(h_BS) - a_hMS;
        d_km_DL = 10^((MAPL_DL - PL_const_DL) / s);
        R_DL(i,j) = d_km_DL * 1000;
        
        % Расчет радиуса UL
        PL_const_UL = A + B*log10(f_ghz*1000) - 13.82*log10(h_BS) - a_hMS;
        d_km_UL = 10^((MAPL_UL - PL_const_UL) / s);
        R_UL(i,j) = d_km_UL * 1000;
    end
end


figure('Position', [200, 200, 1000, 700]);
surf( H, F, R_DL', 'FaceColor', 'red')
hold on;
surf( H,F, R_UL', 'FaceColor', 'blue');
ylabel('Частота, МГц');
xlabel('Высота BS, м');
zlabel('Радиус покрытия, м');
title('Радиусы покрытия Downlink и Uplink');
legend('Downlink', 'Uplink', 'Location', 'northeast');
grid on;
