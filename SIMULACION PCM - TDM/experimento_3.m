%% Sistema TDM para Señales Triangular y Senoidal
% Universidad - Carrera de Telecomunicaciones
% Nivel: Intermedio

clear all; close all; clc;

%% 1. PARÁMETROS DEL SISTEMA
fs = 1000;          % Frecuencia de muestreo [Hz] (debe cumplir Nyquist)
fm_senal = 50;      % Frecuencia máxima de las señales [Hz]
t_final = 0.1;      % Tiempo de simulación [s]
A_triang = 2;       % Amplitud señal triangular [V]
A_sen = 1.5;        % Amplitud señal senoidal [V]

% Verificación del Teorema de Nyquist
if fs <= 2*fm_senal
    error('La frecuencia de muestreo no cumple el teorema de Nyquist');
else
    fprintf('Nyquist verificado: fs = %d Hz > 2*fm = %d Hz\n', fs, 2*fm_senal);
end

%% 2. GENERACIÓN DE SEÑALES ANALÓGICAS ORIGINALES
fprintf('\n=== GENERACIÓN DE SEÑALES ANALÓGICAS ===\n');

% Vector de tiempo de alta resolución para señales analógicas
t_analogico = 0:1/(fs*10):t_final; % 10x mayor resolución que la muestreada

% Señal triangular (diente de sierra)
senal_triangular = A_triang * sawtooth(2*pi*fm_senal*t_analogico, 0.5);
fprintf('Señal triangular: %.1f Hz, %.1f Vpp\n', fm_senal, 2*A_triang);

% Señal senoidal
senal_senoidal = A_sen * sin(2*pi*fm_senal*t_analogico);
fprintf('Señal senoidal: %.1f Hz, %.1f Vpp\n', fm_senal, 2*A_sen);

%% 3. MUESTREO DE LAS SEÑALES (Teorema de Nyquist)
fprintf('\n=== PROCESO DE MUESTREO ===\n');

% Vector de tiempo para muestreo
t_muestreo = 0:1/fs:t_final;
N_muestras = length(t_muestreo);
fprintf('Número de muestras: %d\n', N_muestras);

% Muestreo de las señales
muestras_triangular = A_triang * sawtooth(2*pi*fm_senal*t_muestreo, 0.5);
muestras_senoidal = A_sen * sin(2*pi*fm_senal*t_muestreo);

fprintf('Frecuencia de muestreo: %.0f Hz\n', fs);
fprintf('Periodo de muestreo: %.4f s\n', 1/fs);

%% 4. CUANTIZACIÓN PCM (Pulse Code Modulation)
fprintf('\n=== CUANTIZACIÓN PCM ===\n');

% Parámetros de cuantificación
n_bits = 8;                         % Resolución del cuantificador
niveles = 2^n_bits;                 % Número de niveles de cuantificación
Vmax = max([A_triang, A_sen]) * 1.2; % Rango dinámico del cuantificador
paso_cuant = 2*Vmax / niveles;      % Tamaño del paso de cuantificación

fprintf('Resolución PCM: %d bits\n', n_bits);
fprintf('Niveles de cuantificación: %d\n', niveles);
fprintf('Rango dinámico: ±%.2f V\n', Vmax);
fprintf('Paso de cuantificación: %.4f V\n', paso_cuant);

% Función de cuantificación uniforme
cuantificar = @(x) round(x/paso_cuant) * paso_cuant;

% Señales cuantificadas PCM
pcm_triangular = cuantificar(muestras_triangular);
pcm_senoidal = cuantificar(muestras_senoidal);

% Cálculo del error de cuantificación (ruido)
error_q_triang = mean((muestras_triangular - pcm_triangular).^2);
error_q_sen = mean((muestras_senoidal - pcm_senoidal).^2);

fprintf('Error cuadrático medio cuantificación - Triangular: %.6f V²\n', error_q_triang);
fprintf('Error cuadrático medio cuantificación - Senoidal: %.6f V²\n', error_q_sen);

%% 5. MULTIPLEXACIÓN TDM (Time Division Multiplexing)
fprintf('\n=== MULTIPLEXACIÓN TDM ===\n');

% Creación de la trama TDM
tdm_signal = zeros(1, 2*N_muestras);  % Señal TDM final
tdm_time = 0:1/fs:(2*N_muestras-1)/fs; % Vector de tiempo TDM

% Asignación de slots de tiempo (intercalado)
for i = 1:N_muestras
    tdm_signal(2*i-1) = pcm_triangular(i);  % Slot impar: triangular
    tdm_signal(2*i) = pcm_senoidal(i);      % Slot par: senoidal
end

fprintf('Tasa de símbolos TDM: %.0f baudios\n', 2*fs);
fprintf('Estructura TDM: [Triangular, Senoidal, Triangular, Senoidal, ...]\n');

%% 6. DEMULTIPLEXACIÓN TDM
fprintf('\n=== DEMULTIPLEXACIÓN TDM ===\n');

% Recuperación de las señales desde TDM
recuperado_triangular = zeros(1, N_muestras);
recuperado_senoidal = zeros(1, N_muestras);

for i = 1:N_muestras
    recuperado_triangular(i) = tdm_signal(2*i-1);  % Muestras impares
    recuperado_senoidal(i) = tdm_signal(2*i);      % Muestras pares
end

% Filtrado pasa bajos para reconstrucción (interpolación)
[b, a] = butter(4, 0.8*(fm_senal/(fs/2))); % Filtro Butterworth 4to orden

% Señales reconstruidas
triangular_reconstruida = filter(b, a, recuperado_triangular);
senoidal_reconstruida = filter(b, a, recuperado_senoidal);

fprintf('Señales recuperadas exitosamente\n');

%% 7. ANÁLISIS ESPECTRAL (OPCIONAL - para verificación)
fprintf('\n=== ANÁLISIS ESPECTRAL ===\n');

% Cálculo de FFTs
N_fft = 2^nextpow2(length(t_analogico));
f = fs/2 * linspace(0, 1, N_fft/2+1);

fft_triang_orig = fft(senal_triangular, N_fft);
fft_sen_orig = fft(senal_senoidal, N_fft);
fft_tdm = fft(tdm_signal, N_fft);

%% 8. VISUALIZACIÓN COMPLETA
fprintf('\n=== GENERACIÓN DE GRÁFICAS ===\n');

figure('Position', [100, 100, 1200, 1000]);

% Subplot 1: Señales analógicas originales
subplot(4,2,1);
plot(t_analogico*1000, senal_triangular, 'b-', 'LineWidth', 1.5);
title('Señal Triangular Original (Analógica)');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
grid on;
xlim([0, t_final*1000]);

subplot(4,2,2);
plot(t_analogico*1000, senal_senoidal, 'r-', 'LineWidth', 1.5);
title('Señal Senoidal Original (Analógica)');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
grid on;
xlim([0, t_final*1000]);

% Subplot 2: Señales muestreadas y cuantificadas PCM
subplot(4,2,3);
stem(t_muestreo*1000, muestras_triangular, 'bo', 'MarkerSize', 3);
hold on;
stairs(t_muestreo*1000, pcm_triangular, 'r-', 'LineWidth', 1.5);
plot(t_analogico*1000, senal_triangular, 'k--', 'LineWidth', 0.5);
title('Muestreo y PCM - Señal Triangular');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
legend('Muestras', 'PCM', 'Analógica', 'Location', 'best');
grid on;
xlim([0, t_final*1000]);

subplot(4,2,4);
stem(t_muestreo*1000, muestras_senoidal, 'ro', 'MarkerSize', 3);
hold on;
stairs(t_muestreo*1000, pcm_senoidal, 'b-', 'LineWidth', 1.5);
plot(t_analogico*1000, senal_senoidal, 'k--', 'LineWidth', 0.5);
title('Muestreo y PCM - Señal Senoidal');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
legend('Muestras', 'PCM', 'Analógica', 'Location', 'best');
grid on;
xlim([0, t_final*1000]);

% Subplot 3: Señal TDM
subplot(4,2,5);
stairs(tdm_time*1000, tdm_signal, 'g', 'MarkerSize', 2);
title('Señal Multiplexada TDM');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
grid on;
xlim([0, 20*1000/fs]); % Mostrar primeros 20 periodos

% Subplot 4: Detalle de estructura TDM (primeras muestras)
subplot(4,2,6);
n_detalle = min(10, length(tdm_time));
stairs(tdm_time(1:n_detalle)*1000, tdm_signal(1:n_detalle), 'm', 'MarkerSize', 4);
title('Detalle Estructura TDM (Primeras Muestras)');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
grid on;
for i = 1:n_detalle
    if mod(i,2) == 1
        text(tdm_time(i)*1000, tdm_signal(i)+0.1, 'T', 'HorizontalAlignment', 'center');
    else
        text(tdm_time(i)*1000, tdm_signal(i)+0.1, 'S', 'HorizontalAlignment', 'center');
    end
end

% Subplot 5: Señales recuperadas
subplot(4,2,7);
plot(t_muestreo*1000, recuperado_triangular, 'b-', 'LineWidth', 2);
hold on;
plot(t_analogico*1000, senal_triangular, 'k--', 'LineWidth', 1);
title('Señal Triangular Recuperada vs Original');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
legend('Recuperada', 'Original', 'Location', 'best');
grid on;
xlim([0, t_final*1000]);

subplot(4,2,8);
plot(t_muestreo*1000, recuperado_senoidal, 'r-', 'LineWidth', 2);
hold on;
plot(t_analogico*1000, senal_senoidal, 'k--', 'LineWidth', 1);
title('Señal Senoidal Recuperada vs Original');
xlabel('Tiempo [ms]');
ylabel('Amplitud [V]');
legend('Recuperada', 'Original', 'Location', 'best');
grid on;
xlim([0, t_final*1000]);

sgtitle('Sistema TDM Completo: Generación, PCM, Multiplexación y Recuperación');

%% 9. CÁLCULO DEL ERROR DE RECONSTRUCCIÓN
error_rec_triang = mean((senal_triangular(1:length(t_muestreo)) - recuperado_triangular).^2);
error_rec_sen = mean((senal_senoidal(1:length(t_muestreo)) - recuperado_senoidal).^2);

fprintf('\n=== ERROR DE RECONSTRUCCIÓN ===\n');
fprintf('Error cuadrático medio reconstrucción - Triangular: %.6f V²\n', error_rec_triang);
fprintf('Error cuadrático medio reconstrucción - Senoidal: %.6f V²\n', error_rec_sen);

%% 10. ESPECIFICACIONES TÉCNICAS DEL SISTEMA
fprintf('\n=== ESPECIFICACIONES TÉCNICAS ===\n');
fprintf('Ancho de banda requerido: %.0f Hz\n', fs);
fprintf('Eficiencia espectral TDM: 2 señales/%d Hz = %.3f señales/Hz\n', fs, 2/fs);
fprintf('Retardo total del sistema: %.3f ms\n', 2/fs*1000);

fprintf('\nSimulación completada exitosamente!\n');