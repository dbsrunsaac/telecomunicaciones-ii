% TDM de señales senoidal y triangular con PCM y recuperación
% Autor: Davis (Estudiante de Telecomunicaciones)
% Fecha: 23/09/2025

clc; clear; close all;

%% 1. Generación de señales senoidal y triangular
fs = 1000;              % Frecuencia de muestreo (Hz)
t = 0:1/fs:1;           % Vector de tiempo de 1 segundo

f1 = 5;                 % Frecuencia de la señal senoidal (Hz)
f2 = 5;                 % Frecuencia de la señal triangular (Hz)

senoidal = sin(2*pi*f1*t);              % Señal senoidal
triangular = sawtooth(2*pi*f2*t, 0.5);  % Señal triangular (forma simétrica)

%% 2. Conversión a PCM (Pulse Code Modulation)
n_bits = 8;                             % Número de bits para cuantización
levels = 2^n_bits;                      % Niveles de cuantización

% Normalización entre 0 y 1
senoidal_norm = (senoidal + 1)/2;
triangular_norm = (triangular + 1)/2;

% Cuantización
senoidal_pcm = round(senoidal_norm * (levels - 1));
triangular_pcm = round(triangular_norm * (levels - 1));

%% 3. Muestreo según teorema de Nyquist
% Frecuencia mínima de muestreo: al menos el doble de la frecuencia máxima
fs_nyquist = 2 * max(f1, f2);   % Teorema de Nyquist
ts_nyquist = 1/fs_nyquist;

% Re-muestreo de las señales
t_sampled = 0:ts_nyquist:1;
senoidal_sampled = sin(2*pi*f1*t_sampled);
triangular_sampled = sawtooth(2*pi*f2*t_sampled, 0.5);

% Cuantización de las muestras
senoidal_pcm_sampled = round((senoidal_sampled + 1)/2 * (levels - 1));
triangular_pcm_sampled = round((triangular_sampled + 1)/2 * (levels - 1));

%% 4. Multiplexación TDM (Time Division Multiplexing)
% Intercalamos las muestras de ambas señales
tdm_signal = zeros(1, 2 * length(t_sampled));
tdm_signal(1:2:end) = senoidal_pcm_sampled;
tdm_signal(2:2:end) = triangular_pcm_sampled;

%% 5. Recuperación de señales desde TDM
% Separación de muestras
rec_senoidal = tdm_signal(1:2:end);
rec_triangular = tdm_signal(2:2:end);

% Descuantización (reconstrucción aproximada)
rec_senoidal = (rec_senoidal / (levels - 1)) * 2 - 1;
rec_triangular = (rec_triangular / (levels - 1)) * 2 - 1;

%% Visualización de resultados
figure('Name','TDM y recuperación de señales','NumberTitle','off');

subplot(3,2,1);
plot(t, senoidal); title('Señal Senoidal Original');
xlabel('Tiempo (s)'); ylabel('Amplitud');

subplot(3,2,2);
plot(t, triangular); title('Señal Triangular Original');
xlabel('Tiempo (s)'); ylabel('Amplitud');

subplot(3,2,3);
stem(t_sampled, senoidal_pcm_sampled); title('Senoidal PCM Muestreada');
xlabel('Tiempo (s)'); ylabel('Nivel PCM');

subplot(3,2,4);
stem(t_sampled, triangular_pcm_sampled); title('Triangular PCM Muestreada');
xlabel('Tiempo (s)'); ylabel('Nivel PCM');

subplot(3,2,5);
stem(1:length(tdm_signal), tdm_signal); title('Señal TDM');
xlabel('Índice'); ylabel('Nivel PCM');

subplot(3,2,6);
plot(t_sampled, rec_senoidal, 'r', t_sampled, rec_triangular, 'b');
title('Recuperación de Señales desde TDM');
xlabel('Tiempo (s)'); ylabel('Amplitud');
legend('Senoidal','Triangular');

%% Recomendación adicional
% Para una reconstrucción más precisa, se puede aplicar interpolación o filtrado
% de reconstrucción (por ejemplo, interpolación sinc o filtro pasa bajos).

