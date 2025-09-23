%% Generación de señales triangular y senoidal, PCM, TDM y recuperación
% Este script está diseñado para estudiantes de telecomunicaciones.
% Se generan señales, se convierten a PCM, se multiplexan por división de tiempo (TDM),
% y finalmente se recuperan y visualizan.

clc; % Limpia la ventana de comandos
clear all; % Limpia todas las variables del espacio de trabajo
close all; % Cierra todas las figuras abiertas

%% 1. Generación de señales triangular y senoidal
% Parámetros comunes
fs = 1000; % Frecuencia de muestreo (Hz)
t = 0:1/fs:1; % Vector de tiempo de 0 a 1 segundo

% Señal triangular
f_triangular = 5; % Frecuencia de la señal triangular (Hz)
signal_triangular = sawtooth(2*pi*f_triangular*t, 0.5); % Genera señal triangular

% Señal senoidal
f_senoidal = 10; % Frecuencia de la señal senoidal (Hz)
signal_senoidal = sin(2*pi*f_senoidal*t); % Genera señal senoidal

% Visualización de las señales originales
figure;
subplot(2,1,1);
plot(t, signal_triangular);
title('Señal Triangular Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
plot(t, signal_senoidal);
title('Señal Senoidal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

%% 2. Conversión de las señales a PCM
% Parámetros para PCM
bits = 8; % Número de bits para la cuantización
L = 2^bits; % Número de niveles de cuantización
Vmax = max([max(signal_triangular), max(signal_senoidal)]); % Valor máximo de las señales
Vmin = min([min(signal_triangular), min(signal_senoidal)]); % Valor mínimo de las señales
delta = (Vmax - Vmin) / L; % Paso de cuantización

% Cuantización de la señal triangular
signal_triangular_pcm = floor((signal_triangular - Vmin) / delta) * delta + Vmin + delta/2;

% Cuantización de la señal senoidal
signal_senoidal_pcm = floor((signal_senoidal - Vmin) / delta) * delta + Vmin + delta/2;

% Visualización de las señales cuantizadas (PCM)
figure;
subplot(2,1,1);
stairs(t, signal_triangular_pcm);
title('Señal Triangular Cuantizada (PCM)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
stairs(t, signal_senoidal_pcm);
title('Señal Senoidal Cuantizada (PCM)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

%% 3. Emulación de un muestreador respetando el teorema de Nyquist
% Frecuencia de muestreo para TDM (debe ser al menos el doble de la frecuencia máxima de las señales)
fs_muestreo = 2 * max(f_triangular, f_senoidal) * 10; % Se multiplica por 10 para evitar aliasing
t_muestreo = 0:1/fs_muestreo:1; % Vector de tiempo para el muestreo

% Muestreo de las señales
signal_triangular_muestreada = interp1(t, signal_triangular_pcm, t_muestreo);
signal_senoidal_muestreada = interp1(t, signal_senoidal_pcm, t_muestreo);

% Visualización de las señales muestreadas
figure;
subplot(2,1,1);
stem(t_muestreo, signal_triangular_muestreada, 'filled');
title('Señal Triangular Muestreada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
stem(t_muestreo, signal_senoidal_muestreada, 'filled');
title('Señal Senoidal Muestreada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

%% 4. Obtención de la señal TDM
% Intercalación de las muestras de ambas señales para formar la señal TDM
signal_tdm = zeros(1, 2*length(t_muestreo));
signal_tdm(1:2:end) = signal_triangular_muestreada; % Muestras de la señal triangular en posiciones impares
signal_tdm(2:2:end) = signal_senoidal_muestreada; % Muestras de la señal senoidal en posiciones pares

% Visualización de la señal TDM
figure;
stem((1:length(signal_tdm))/fs_muestreo, signal_tdm, 'filled');
title('Señal TDM');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

%% 5. Recuperación de las señales transmitidas por TDM
% Separación de las muestras en la señal TDM
signal_triangular_recuperada = signal_tdm(1:2:end); % Recupera la señal triangular
signal_senoidal_recuperada = signal_tdm(2:2:end); % Recupera la señal senoidal

% Reconstrucción de las señales recuperadas
t_recuperado = 0:1/fs_muestreo:1-1/fs_muestreo;
signal_triangular_recuperada = interp1(t_recuperado, signal_triangular_recuperada, t, 'linear');
signal_senoidal_recuperada = interp1(t_recuperado, signal_senoidal_recuperada, t, 'linear');

% Visualización de las señales recuperadas
figure;
subplot(2,1,1);
plot(t, signal_triangular_recuperada);
title('Señal Triangular Recuperada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
plot(t, signal_senoidal_recuperada);
title('Señal Senoidal Recuperada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

%% Visualización final de todas las señales en una figura
figure;
subplot(4,1,1);
plot(t, signal_triangular);
title('Señal Triangular Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(4,1,2);
plot(t, signal_triangular_recuperada);
title('Señal Triangular Recuperada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(4,1,3);
plot(t, signal_senoidal);
title('Señal Senoidal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(4,1,4);
plot(t, signal_senoidal_recuperada);
title('Señal Senoidal Recuperada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
