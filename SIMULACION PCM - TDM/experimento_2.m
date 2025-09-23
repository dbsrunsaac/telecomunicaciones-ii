%% SISTEMA DE MODULACIÓN PCM Y MULTIPLEXACIÓN TDM
% Autor: Sistema de comunicaciones digitales
% Descripción: Este código implementa la conversión PCM y multiplexación TDM
% de señales triangular y senoidal, con recuperación de las señales originales.
% Nivel: Intermedio - Estudiantes de telecomunicaciones

clear all; close all; clc;

%% 1. PARÁMETROS INICIALES Y GENERACIÓN DE SEÑALES
% Definimos los parámetros temporales y de frecuencia
fs_original = 10000;    % Frecuencia de muestreo original (Hz) - alta para simular señal "continua"
duracion = 0.01;        % Duración de las señales (segundos)
t = 0:1/fs_original:duracion; % Vector de tiempo

% Frecuencias de las señales (deben cumplir Nyquist)
f_triangular = 100;     % Frecuencia de la señal triangular (Hz)
f_senoidal = 150;       % Frecuencia de la señal senoidal (Hz)

% GENERACIÓN DE SEÑAL TRIANGULAR
% sawtooth() genera una onda diente de sierra, con parámetro 0.5 se vuelve triangular
% El factor 2*pi*f convierte la frecuencia a radianes/segundo
senal_triangular = sawtooth(2*pi*f_triangular*t, 0.5);

% GENERACIÓN DE SEÑAL SENOIDAL
% sin() genera una señal senoidal pura
senal_senoidal = sin(2*pi*f_senoidal*t);

% Normalización de amplitudes para PCM
% Esto asegura que las señales estén en el rango [-1, 1]
senal_triangular = senal_triangular / max(abs(senal_triangular));
senal_senoidal = senal_senoidal / max(abs(senal_senoidal));

%% 2. CONVERSIÓN A PCM (PULSE CODE MODULATION)
% PCM convierte señales analógicas a digitales mediante:
% - Muestreo
% - Cuantización
% - Codificación

% Parámetros PCM
n_bits = 8;                     % Número de bits para cuantización
L = 2^n_bits;                   % Número de niveles de cuantización
Vmax = 1;                       % Voltaje máximo
Vmin = -1;                      % Voltaje mínimo
delta = (Vmax - Vmin) / L;      % Paso de cuantización

% Función para realizar PCM
function [pcm_signal, niveles_cuant] = realizar_pcm(signal, L, Vmin, Vmax)
    % CUANTIZACIÓN: Asigna cada muestra al nivel más cercano
    delta = (Vmax - Vmin) / L;
    niveles_cuant = Vmin:delta:Vmax-delta;
    
    % Inicializar señal PCM
    pcm_signal = zeros(size(signal));
    
    % Cuantizar cada muestra
    for i = 1:length(signal)
        % Encuentra el nivel de cuantización más cercano
        [~, indice] = min(abs(signal(i) - niveles_cuant));
        pcm_signal(i) = niveles_cuant(indice);
    end
end

% Aplicar PCM a ambas señales
[pcm_triangular, niveles] = realizar_pcm(senal_triangular, L, Vmin, Vmax);
[pcm_senoidal, ~] = realizar_pcm(senal_senoidal, L, Vmin, Vmax);

%% 3. MUESTREO SEGÚN TEOREMA DE NYQUIST
% El teorema de Nyquist establece que fs >= 2*fmax para evitar aliasing

% Determinar frecuencia de muestreo mínima
f_max = max(f_triangular, f_senoidal); % Frecuencia máxima presente
fs_nyquist = 2.5 * f_max;              % Factor de seguridad de 2.5x
fprintf('Frecuencia máxima de señal: %.2f Hz\n', f_max);
fprintf('Frecuencia de muestreo Nyquist (mínima): %.2f Hz\n', 2*f_max);
fprintf('Frecuencia de muestreo utilizada: %.2f Hz\n', fs_nyquist);

% Período de muestreo
Ts = 1/fs_nyquist;

% EMULACIÓN DEL MUESTREADOR
% Tomamos muestras a intervalos regulares Ts
indices_muestreo = 1:round(fs_original/fs_nyquist):length(t);
t_muestreado = t(indices_muestreo);

% Muestras de las señales PCM
muestras_pcm_triangular = pcm_triangular(indices_muestreo);
muestras_pcm_senoidal = pcm_senoidal(indices_muestreo);

%% 4. MULTIPLEXACIÓN TDM (TIME DIVISION MULTIPLEXING)
% TDM intercala muestras de diferentes señales en el tiempo

% Crear señal TDM intercalando muestras
% Método: [muestra1_triangular, muestra1_senoidal, muestra2_triangular, muestra2_senoidal, ...]
num_muestras = length(muestras_pcm_triangular);
señal_tdm = zeros(1, 2*num_muestras);

% Intercalar las muestras (multiplexación)
señal_tdm(1:2:end) = muestras_pcm_triangular;  % Posiciones impares: triangular
señal_tdm(2:2:end) = muestras_pcm_senoidal;    % Posiciones pares: senoidal

% Vector de tiempo para TDM
% El período de símbolo en TDM es la mitad del período de muestreo original
Ts_tdm = Ts / 2;
t_tdm = 0:Ts_tdm:(length(señal_tdm)-1)*Ts_tdm;

%% 5. DEMULTIPLEXACIÓN Y RECUPERACIÓN DE SEÑALES
% Separar las señales del flujo TDM

% DEMULTIPLEXACIÓN: Extraer muestras alternadas
recuperada_triangular_muestras = señal_tdm(1:2:end);  % Muestras impares
recuperada_senoidal_muestras = señal_tdm(2:2:end);    % Muestras pares

% INTERPOLACIÓN para reconstruir señales continuas
% Usamos interpolación spline para una reconstrucción suave
t_recuperado = 0:1/fs_original:t_muestreado(end);

% Interpolación de señales recuperadas
if length(t_muestreado) > 1 && length(recuperada_triangular_muestras) > 1
    recuperada_triangular = interp1(t_muestreado, recuperada_triangular_muestras, ...
                                   t_recuperado, 'spline');
    recuperada_senoidal = interp1(t_muestreado, recuperada_senoidal_muestras, ...
                                 t_recuperado, 'spline');
else
    recuperada_triangular = recuperada_triangular_muestras;
    recuperada_senoidal = recuperada_senoidal_muestras;
end

% Aplicar filtro paso bajo para suavizar (simula filtro de reconstrucción)
% butterworth es un filtro IIR con respuesta plana en banda de paso
fc_filtro = f_max * 1.2;  % Frecuencia de corte del filtro
[b, a] = butter(4, fc_filtro/(fs_original/2), 'low');  % Filtro Butterworth orden 4
recuperada_triangular_filtrada = filtfilt(b, a, recuperada_triangular);
recuperada_senoidal_filtrada = filtfilt(b, a, recuperada_senoidal);

%% 6. VISUALIZACIÓN DE RESULTADOS
figure('Position', [100, 100, 1400, 900]);
sgtitle('Sistema PCM-TDM: Transmisión y Recuperación de Señales', 'FontSize', 14, 'FontWeight', 'bold');

% Subplot 1: Señal triangular original y PCM
subplot(3,3,1);
plot(t*1000, senal_triangular, 'b-', 'LineWidth', 1.5);
hold on;
plot(t*1000, pcm_triangular, 'r--', 'LineWidth', 1);
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title('Señal Triangular Original vs PCM');
legend('Original', 'PCM', 'Location', 'best');
grid on;

% Subplot 2: Señal senoidal original y PCM
subplot(3,3,2);
plot(t*1000, senal_senoidal, 'b-', 'LineWidth', 1.5);
hold on;
plot(t*1000, pcm_senoidal, 'r--', 'LineWidth', 1);
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title('Señal Senoidal Original vs PCM');
legend('Original', 'PCM', 'Location', 'best');
grid on;

% Subplot 3: Espectro de frecuencias
subplot(3,3,3);
N_fft = 2048;
f_axis = linspace(0, fs_original/2, N_fft/2);
fft_triangular = abs(fft(senal_triangular, N_fft));
fft_senoidal = abs(fft(senal_senoidal, N_fft));
plot(f_axis, fft_triangular(1:N_fft/2), 'b-', 'LineWidth', 1);
hold on;
plot(f_axis, fft_senoidal(1:N_fft/2), 'r-', 'LineWidth', 1);
xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
title('Espectros de Frecuencia');
legend('Triangular', 'Senoidal');
xlim([0 500]);
grid on;

% Subplot 4: Muestras según Nyquist - Triangular
subplot(3,3,4);
stem(t_muestreado*1000, muestras_pcm_triangular, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title(['Muestras PCM Triangular (fs = ' num2str(fs_nyquist, '%.1f') ' Hz)']);
grid on;

% Subplot 5: Muestras según Nyquist - Senoidal
subplot(3,3,5);
stem(t_muestreado*1000, muestras_pcm_senoidal, 'r', 'LineWidth', 1.5);
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title(['Muestras PCM Senoidal (fs = ' num2str(fs_nyquist, '%.1f') ' Hz)']);
grid on;

% Subplot 6: Señal TDM multiplexada
subplot(3,3,6);
stairs(t_tdm*1000, señal_tdm, 'k', 'LineWidth', 1.5);
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title('Señal TDM Multiplexada');
grid on;
% Añadir marcadores para identificar origen de muestras
hold on;
for i = 1:min(20, length(señal_tdm))  % Mostrar solo primeras 20 muestras
    if mod(i, 2) == 1
        plot(t_tdm(i)*1000, señal_tdm(i), 'bo', 'MarkerSize', 6);
    else
        plot(t_tdm(i)*1000, señal_tdm(i), 'ro', 'MarkerSize', 6);
    end
end
if length(señal_tdm) <= 20
    legend('TDM', 'Triangular', 'Senoidal', 'Location', 'best');
end

% Subplot 7: Señal triangular recuperada
subplot(3,3,7);
plot(t*1000, senal_triangular, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original');
hold on;
plot(t_recuperado*1000, recuperada_triangular_filtrada, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Recuperada');
stem(t_muestreado*1000, recuperada_triangular_muestras, 'r', 'LineWidth', 0.5, 'MarkerSize', 4, 'DisplayName', 'Muestras');
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title('Recuperación Señal Triangular');
legend('Location', 'best');
grid on;

% Subplot 8: Señal senoidal recuperada
subplot(3,3,8);
plot(t*1000, senal_senoidal, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original');
hold on;
plot(t_recuperado*1000, recuperada_senoidal_filtrada, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Recuperada');
stem(t_muestreado*1000, recuperada_senoidal_muestras, 'r', 'LineWidth', 0.5, 'MarkerSize', 4, 'DisplayName', 'Muestras');
xlabel('Tiempo (ms)'); ylabel('Amplitud');
title('Recuperación Señal Senoidal');
legend('Location', 'best');
grid on;

% Subplot 9: Métricas de calidad
subplot(3,3,9);
% Calcular SNR (Signal-to-Noise Ratio)
% Alinear señales para comparación
min_length = min(length(senal_triangular), length(recuperada_triangular_filtrada));
error_triangular = senal_triangular(1:min_length) - recuperada_triangular_filtrada(1:min_length);
error_senoidal = senal_senoidal(1:min_length) - recuperada_senoidal_filtrada(1:min_length);

SNR_triangular = 10*log10(sum(senal_triangular(1:min_length).^2) / sum(error_triangular.^2));
SNR_senoidal = 10*log10(sum(senal_senoidal(1:min_length).^2) / sum(error_senoidal.^2));

% Mostrar métricas
text(0.1, 0.9, 'MÉTRICAS DEL SISTEMA', 'FontWeight', 'bold', 'FontSize', 12);
text(0.1, 0.7, sprintf('Bits PCM: %d', n_bits));
text(0.1, 0.6, sprintf('Niveles cuantización: %d', L));
text(0.1, 0.5, sprintf('Paso cuantización: %.4f', delta));
text(0.1, 0.4, sprintf('Frecuencia muestreo: %.1f Hz', fs_nyquist));
text(0.1, 0.3, sprintf('SNR Triangular: %.2f dB', SNR_triangular));
text(0.1, 0.2, sprintf('SNR Senoidal: %.2f dB', SNR_senoidal));
axis off;

%% ANÁLISIS ADICIONAL Y COMENTARIOS EDUCATIVOS

fprintf('\n=== ANÁLISIS DEL SISTEMA PCM-TDM ===\n');
fprintf('1. TEOREMA DE NYQUIST:\n');
fprintf('   - Para evitar aliasing: fs >= 2*fmax\n');
fprintf('   - Frecuencia máxima en el sistema: %.2f Hz\n', f_max);
fprintf('   - Frecuencia de muestreo mínima: %.2f Hz\n', 2*f_max);
fprintf('   - Frecuencia utilizada (con margen): %.2f Hz\n\n', fs_nyquist);

fprintf('2. MODULACIÓN PCM:\n');
fprintf('   - Resolución: %d bits\n', n_bits);
fprintf('   - Niveles de cuantización: %d\n', L);
fprintf('   - Rango dinámico: %.2f dB\n', 20*log10(L));
fprintf('   - Error de cuantización máximo: ±%.4f\n\n', delta/2);

fprintf('3. MULTIPLEXACIÓN TDM:\n');
fprintf('   - Número de canales: 2\n');
fprintf('   - Velocidad de símbolos TDM: %.2f símbolos/s\n', 1/Ts_tdm);
fprintf('   - Ancho de banda requerido: %.2f Hz (aproximado)\n\n', 1/(2*Ts_tdm));

fprintf('4. CALIDAD DE RECUPERACIÓN:\n');
fprintf('   - SNR Señal Triangular: %.2f dB\n', SNR_triangular);
fprintf('   - SNR Señal Senoidal: %.2f dB\n', SNR_senoidal);
fprintf('   - Distorsión principalmente por: cuantización PCM y filtrado\n');

%% ESPECIFICACIONES ADICIONALES PARA PERSONALIZACIÓN

fprintf('\n=== ESPECIFICACIONES PARA MODIFICAR EL CÓDIGO ===\n');
fprintf('Para adaptar este código a diferentes requisitos, modifique:\n\n');
fprintf('1. SEÑALES DE ENTRADA:\n');
fprintf('   - f_triangular, f_senoidal: Frecuencias de las señales (líneas 14-15)\n');
fprintf('   - duracion: Duración de la simulación (línea 10)\n\n');
fprintf('2. PARÁMETROS PCM:\n');
fprintf('   - n_bits: Resolución de cuantización (línea 35)\n');
fprintf('   - Vmax, Vmin: Rango de voltajes (líneas 37-38)\n\n');
fprintf('3. MUESTREO:\n');
fprintf('   - Factor en fs_nyquist: Margen sobre Nyquist (línea 71)\n\n');
fprintf('4. FILTRADO:\n');
fprintf('   - Orden del filtro Butterworth (línea 128)\n');
fprintf('   - Factor de frecuencia de corte (línea 127)\n\n');
fprintf('5. VISUALIZACIÓN:\n');
fprintf('   - N_fft: Resolución del espectro (línea 162)\n');
fprintf('   - Número de subplots y su disposición\n\n');
fprintf('6. PARA AÑADIR MÁS SEÑALES:\n');
fprintf('   - Crear nueva señal en sección 1\n');
fprintf('   - Aplicar PCM en sección 2\n');
fprintf('   - Muestrear en sección 3\n');
fprintf('   - Modificar intercalado TDM en sección 4 (usar mod para N señales)\n');
fprintf('   - Ajustar demultiplexación en sección 5\n');