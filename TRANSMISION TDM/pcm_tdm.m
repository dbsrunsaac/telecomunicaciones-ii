%% TRANSMISIÓN DE SEÑALES MEDIANTE PCM Y TDM
% Autor: Asistente de Telecomunicaciones
% Fecha: 2023
% Descripción: Emulación de transmisión PCM-TDM para señales triangular y senoidal

clc; clear all; close all;

%% PARÁMETROS DEL SISTEMA
% =========================================================================
% Estos parámetros definen las características del sistema de transmisión

% Parámetros de las señales
f_senoidal = 100;          % Frecuencia de la señal senoidal [Hz]
f_triangular = 50;         % Frecuencia de la señal triangular [Hz]
Amplitud = 1;              % Amplitud máxima de las señales [V]
duracion = 0.02;           % Duración de la simulación [s]

% Parámetros de muestreo (Nyquist)
f_muestreo = 2000;         % Frecuencia de muestreo [Hz]
t_muestreo = 1/f_muestreo; % Período de muestreo [s]

% Parámetros PCM
niveles_cuantizacion = 8;  % Número de niveles de cuantización (3 bits)
bits_por_muestra = 3;      % Bits por muestra

% Parámetros TDM
num_canales = 2;           % Número de canales a multiplexar

fprintf('=== SISTEMA DE TRANSMISIÓN PCM-TDM ===\n');
fprintf('Frecuencia de muestreo: %d Hz\n', f_muestreo);
fprintf('Bits por muestra: %d bits\n', bits_por_muestra);
fprintf('Niveles de cuantización: %d niveles\n', niveles_cuantizacion);
fprintf('Número de canales: %d\n\n', num_canales);

%% GENERACIÓN DE SEÑALES ANALÓGICAS
% =========================================================================
% Creación del vector de tiempo
t = 0:t_muestreo:duracion;

% Señal senoidal
senal_senoidal = Amplitud * sin(2*pi*f_senoidal*t);

% Señal triangular (usando la función sawtooth)
senal_triangular = Amplitud * sawtooth(2*pi*f_triangular*t + pi/2, 0.5);

% Visualización de las señales originales
figure('Name', 'Señales Analógicas Originales', 'Position', [100 100 900 600]);

subplot(2,1,1);
plot(t, senal_senoidal, 'b-', 'LineWidth', 1.5);
title('Señal Senoidal Original');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
grid on;
xlim([0 duracion]);

subplot(2,1,2);
plot(t, senal_triangular, 'r-', 'LineWidth', 1.5);
title('Señal Triangular Original');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
grid on;
xlim([0 duracion]);

%% PROCESO DE MUESTREO
% =========================================================================
% El muestreo convierte la señal continua en discreta
% Se utiliza la misma frecuencia de muestreo para ambas señales

fprintf('Proceso de muestreo completado:\n');
fprintf('- Número de muestras: %d\n', length(t));
fprintf('- Frecuencia de muestreo: %.2f Hz\n', f_muestreo);
fprintf('- Período de muestreo: %.6f s\n\n', t_muestreo);

%% PROCESO DE CUANTIZACIÓN PCM
% =========================================================================
% La cuantización asigna valores discretos a las muestras analógicas

% Crear niveles de cuantización uniforme
niveles = linspace(-Amplitud, Amplitud, niveles_cuantizacion);

% Función para cuantizar una señal
cuantizar = @(senal) interp1(niveles, niveles, senal, 'nearest', 'extrap');

% Cuantización de las señales
senal_senoidal_cuant = cuantizar(senal_senoidal);
senal_triangular_cuant = cuantizar(senal_triangular);

% Visualización de la cuantización
figure('Name', 'Proceso de Cuantización PCM', 'Position', [200 100 900 600]);

subplot(2,1,1);
stairs(t, senal_senoidal_cuant, 'b-', 'LineWidth', 1.5); hold on;
plot(t, senal_senoidal, 'r--', 'LineWidth', 1);
title('Cuantización PCM - Señal Senoidal');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
legend('Señal Cuantizada', 'Señal Original');
grid on;
xlim([0 0.005]);

subplot(2,1,2);
stairs(t, senal_triangular_cuant, 'b-', 'LineWidth', 1.5); hold on;
plot(t, senal_triangular, 'r--', 'LineWidth', 1);
title('Cuantización PCM - Señal Triangular');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
legend('Señal Cuantizada', 'Señal Original');
grid on;
xlim([0 0.005]);

%% CODIFICACIÓN PCM
% =========================================================================
% Convertir los niveles cuantizados a código binario

% Función para codificar a binario
codificar_binario = @(senal) de2bi(round((senal + Amplitud) * (niveles_cuantizacion-1) / (2*Amplitud)), bits_por_muestra, 'left-msb');

% Codificar las señales
codigo_senoidal = codificar_binario(senal_senoidal_cuant);
codigo_triangular = codificar_binario(senal_triangular_cuant);

fprintf('Codificación PCM completada:\n');
fprintf('- Tasa de bits: %d bps\n', f_muestreo * bits_por_muestra);
fprintf('- Muestra 1 Senoidal: %s\n', num2str(codigo_senoidal(1,:)));
fprintf('- Muestra 1 Triangular: %s\n\n', num2str(codigo_triangular(1,:)));

%% MULTIPLEXACIÓN POR DIVISIÓN DE TIEMPO (TDM)
% =========================================================================
% Intercalar las muestras de ambos canales en el tiempo

% Reorganizar los códigos binarios para TDM
num_muestras = size(codigo_senoidal, 1);
senal_tdm = zeros(num_muestras * num_canales, bits_por_muestra);

for i = 1:num_muestras
    % Canal 1 (Senoidal)
    senal_tdm(2*i-1, :) = codigo_senoidal(i, :);
    
    % Canal 2 (Triangular)
    senal_tdm(2*i, :) = codigo_triangular(i, :);
end

fprintf('Multiplexación TDM completada:\n');
fprintf('- Trama TDM completa: %d bits\n', numel(senal_tdm));
fprintf('- Primera trama TDM: %s%s\n\n', ...
    num2str(senal_tdm(1,:)), num2str(senal_tdm(2,:)));

%% VISUALIZACIÓN DE LA SEÑAL TDM
% =========================================================================
% Crear vector de tiempo para TDM
t_tdm = (0:size(senal_tdm,1)-1) * (t_muestreo/num_canales);

figure('Name', 'Señal Multiplexada TDM', 'Position', [300 100 900 400]);

% Convertir matriz binaria a señal digital
senal_digital = zeros(size(senal_tdm,1), 1);
for i = 1:size(senal_tdm,1)
    % Convertir binario a decimal para visualización
    senal_digital(i) = bi2de(senal_tdm(i,:), 'left-msb');
end

stairs(t_tdm, senal_digital, 'b-', 'LineWidth', 1.5);
title('Señal Multiplexada TDM (Valores Decimales de las Palabras PCM)');
xlabel('Tiempo [s]');
ylabel('Valor de Palabra PCM');
grid on;
xlim([0 0.005]);

%% DEMULTIPLEXACIÓN TDM
% =========================================================================
% Separar la señal TDM en los canales originales

senal_demux_senoidal = zeros(num_muestras, bits_por_muestra);
senal_demux_triangular = zeros(num_muestras, bits_por_muestra);

for i = 1:num_muestras
    % Extraer canal 1 (Senoidal)
    senal_demux_senoidal(i, :) = senal_tdm(2*i-1, :);
    
    % Extraer canal 2 (Triangular)
    senal_demux_triangular(i, :) = senal_tdm(2*i, :);
end

fprintf('Demultiplexación TDM completada:\n');
fprintf('- Canales recuperados correctamente\n\n');

%% DECODIFICACIÓN PCM
% =========================================================================
% Convertir binario a niveles de cuantización

% Función para decodificar
decodificar_binario = @(codigo) (bi2de(codigo, 'left-msb') * (2*Amplitud) / (niveles_cuantizacion-1)) - Amplitud;

% Decodificar las señales
senal_recuperada_senoidal = decodificar_binario(senal_demux_senoidal);
senal_recuperada_triangular = decodificar_binario(senal_demux_triangular);

%% ANÁLISIS DE CALIDAD
% =========================================================================
% Calcular error cuadrático medio (MSE)

mse_senoidal = mean((senal_senoidal_cuant - senal_recuperada_senoidal).^2);
mse_triangular = mean((senal_triangular_cuant - senal_recuperada_triangular).^2);

fprintf('Análisis de Calidad:\n');
fprintf('- MSE Senoidal: %.6f\n', mse_senoidal);
fprintf('- MSE Triangular: %.6f\n\n', mse_triangular);

%% VISUALIZACIÓN FINAL
% =========================================================================
figure('Name', 'Comparación: Original vs Recuperada', 'Position', [400 100 900 600]);

subplot(2,1,1);
stairs(t, senal_recuperada_senoidal, 'b-', 'LineWidth', 1.5); hold on;
plot(t, senal_senoidal, 'r--', 'LineWidth', 1);
title('Señal Senoidal: Original vs Recuperada');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
legend('Recuperada (PCM)', 'Original');
grid on;
xlim([0 0.01]);

subplot(2,1,2);
stairs(t, senal_recuperada_triangular, 'b-', 'LineWidth', 1.5); hold on;
plot(t, senal_triangular, 'r--', 'LineWidth', 1);
title('Señal Triangular: Original vs Recuperada');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
legend('Recuperada (PCM)', 'Original');
grid on;
xlim([0 0.01]);

%% ESPECTRO DE LA SEÑAL TDM
% =========================================================================
figure('Name', 'Espectro de la Señal TDM', 'Position', [500 100 900 400]);

% Calcular FFT
L = length(senal_digital);
Y = fft(senal_digital);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = f_muestreo * num_canales * (0:(L/2))/L;

plot(f, P1, 'b-', 'LineWidth', 1.5);
title('Espectro de Frecuencia de la Señal TDM');
xlabel('Frecuencia [Hz]');
ylabel('|P1(f)|');
grid on;
xlim([0 f_muestreo * num_canales / 2]);

fprintf('=== SIMULACIÓN COMPLETADA ===\n');