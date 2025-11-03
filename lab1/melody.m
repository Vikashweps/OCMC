Fs = 44100;
amplitude = 1;
%temp = 120;
%beat_duration = 60/temp;

% Частоты нот 
notes_freq.G3 = 196.00;    % Соль 
notes_freq.A3 = 220.00;    % Ля 
notes_freq.B3 = 246.94;
notes_freq.C4 = 261.63;    % До 
notes_freq.D4 = 293.66;    % Ре 
notes_freq.E4 = 329.63;    % Ми 
notes_freq.F4 = 349.23;    % Фа 
notes_freq.G4 = 392.00;    % Соль 
notes_freq.A4 = 440.00;    % Ля 
notes_freq.B4 = 493.88;
notes_freq.C5 = 523.25;
notes_freq.D5 = 587.33;
melody_sequence = {
    'G3', 1;   % Соль 
    'E4', 1;   % ми 
    'E4', 1;   % ми 
    'D4', 1;   % ре 
    'E4', 1;   % ми 
    'C4', 1;   % до 
    'G3', 1;   % соль 
    'G3', 1;   % соль 
    
    'G3', 1;   % соль 
    'E4', 1;   % ми 
    'E4', 1;   % ми 
    'F4', 1;   % фа 
    'D4', 1;   % ре 
    'G4', 1;   % соль 
    
    'G4', 1;   % Соль 
    'A3', 1;   % ля 
    'A3', 1;   % ля 
    'F4', 1;   % фа 
    'F4', 1;   % фа 
    'E4', 1;   % ми 
    'D4', 1;   % ре 
    'C4', 1;   % до 
    'A3', 1;   % ля 
    'E4', 1;   % ми 
    'E4', 1;   % ми 
    'D4', 1;   % ре 
    'E4', 1;   % ми 
    'C4', 1;   % до 
};

% Генерация мелодии
full_melody = [];
for i = 1:length(melody_sequence)
    note_name = melody_sequence{i, 1};
    %duration_beats = melody_sequence{i, 2};
    note_duration = melody_sequence{i, 2} * 0.5;
    freq = notes_freq.(note_name);
    
    t = 0:1/Fs:note_duration;
    tone = amplitude * sin(2*pi*freq*t);
    
    full_melody = [full_melody, tone];
end

% Нормализация
full_melody = full_melody / max(abs(full_melody));

player = audioplayer(full_melody, Fs);
play(player);
fprintf('Воспроизведение началось. Длительность: %.1f сек\n', length(full_melody)/Fs);
fprintf('Чтобы остановить: stop(player)\n');
