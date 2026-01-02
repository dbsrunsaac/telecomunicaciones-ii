clc; clear; close all;

%% Parámetros de la simulación
order = 32;
symorder = 'gray';

symbols = 0 : order - 1;
bits = reshape(de2bi(symbols).', [], 1);

%% Modulación QAM
qam_symbols = qammod(bits, order, symorder, "InputType", "bit", "UnitAveragePower", true);
n_symbols = length(qam_symbols);

%% Parámetros de ruido
SNR_dB = 20;  % Relación señal-ruido en dB (ajustable)
% Valores típicos:
% 30 dB: muy poco ruido
% 20 dB: ruido moderado
% 10 dB: mucho ruido
% 5 dB: ruido severo

%% Plot diagrama de constelación ORIGINAL
axlim = 1.2;
color_original = '#FF5757';
figure('Position', [100 100 1200 500]);
subplot(1,2,1);
hold on;
for i = 0 : n_symbols - 1
    txt = bits(i * log2(order) + 1 : i * log2(order) + log2(order));
    txt = num2str(txt.');
    txt = strrep(txt, ' ', '');
    text(real(qam_symbols(i + 1)) - 0.1, imag(qam_symbols(i + 1)) + 0.1, txt, 'FontSize', 7);
end
plot(real(qam_symbols), imag(qam_symbols), '.', 'MarkerSize', 30, 'Color', color_original);
hold off;
ylim([-axlim, axlim]); xlim([-axlim axlim]);
set(gca, 'XAxisLocation', 'origin'); set(gca, 'YAxisLocation', 'origin');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
yticks([]); xticks([]); axis square;
title('Constelación 32QAM Original');
grid on;

%% Señal modulada
sampling_frequency = 1000;
time = (0 : n_symbols * sampling_frequency - 1) / sampling_frequency;

carrier_I = cos(2 * pi * time);
carrier_Q = sin(2 * pi * time);

qam_symbols_repeated = repelem(qam_symbols, sampling_frequency).';
modulated_signal = real(qam_symbols_repeated) .* carrier_I + imag(qam_symbols_repeated) .* carrier_Q;

%% AGREGAR RUIDO A LA SEÑAL MODULADA
% Calcular la potencia de la señal
signal_power = mean(abs(modulated_signal).^2);

% Convertir SNR de dB a escala lineal
SNR_linear = 10^(SNR_dB/10);

% Calcular la potencia del ruido
noise_power = signal_power / SNR_linear;

% Generar ruido blanco gaussiano (AWGN)
noise = sqrt(noise_power/2) * (randn(size(modulated_signal)) + 1i*randn(size(modulated_signal)));
noise_real = real(noise);  % Usamos solo la parte real ya que la señal modulada es real

% Señal con ruido
noisy_signal = modulated_signal + noise_real;

%% Demodulación para obtener constelación con ruido
% Para obtener los símbolos con ruido, necesitamos demodular
time_symbols = 0:sampling_frequency:length(modulated_signal)-1;
symbol_indices = floor(time_symbols/sampling_frequency) + 1;

% Extraer símbolos en los instantes de muestreo
received_symbols = qam_symbols_repeated(symbol_indices);

% Agregar ruido a los símbolos (en baseband)
symbol_noise_power = mean(abs(qam_symbols).^2) / SNR_linear;
symbol_noise = sqrt(symbol_noise_power/2) * (randn(size(qam_symbols)) + 1i*randn(size(qam_symbols)));
noisy_symbols = qam_symbols + symbol_noise;

%% Plot diagrama de constelación CON RUIDO
subplot(1,2,2);
hold on;
plot(real(noisy_symbols), imag(noisy_symbols), '.', 'MarkerSize', 20, 'Color', '#233ce6');
plot(real(qam_symbols), imag(qam_symbols), 'r.', 'MarkerSize', 5, 'Color', '#FF5757');
hold off;
ylim([-axlim, axlim]); xlim([-axlim axlim]);
set(gca, 'XAxisLocation', 'origin'); set(gca, 'YAxisLocation', 'origin');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
yticks([]); xticks([]); axis square;
title(sprintf('Constelación 32QAM con Ruido (SNR = %d dB)', SNR_dB));
grid on;
legend('Símbolos con ruido', 'Posiciones ideales', 'Location', 'best');

%% Plot señales moduladas
color = '#233ce6';
figure('Position', [100 100 800 800]);
tiledlayout(4, 1);

% Señal de bits
nexttile;
x_values = 0 : length(bits);
x_values = repelem(x_values, 2);
x_values = x_values(2 : end - 1);
y_values = repelem(bits, 2);
plot(x_values, y_values, 'Color', color, 'LineWidth', 1.5);
max_ticks = 16; step = length(bits) / max_ticks;
if step < 1
    ticks = 0 : length(bits);
else
    ticks = 0 : step : length(bits);
end
xticks(ticks); xlabel('bits');
ylim([-0.2, 1.2]); xlim([0, x_values(length(x_values))]);
title('Señal Digital Original');
grid on;

% Señal modulada sin ruido
nexttile;
plot(time, modulated_signal, 'Color', color, 'LineWidth', 1.5);
max_ticks = 16; step = n_symbols / max_ticks;
if step < 1
    ticks = 0 : n_symbols;
else
    ticks = 0 : step : n_symbols;
end
xticks(ticks); xlim([0, n_symbols]); xlabel('símbolos');
title('Señal Modulada 32QAM (sin ruido)');
grid on;

% Ruido
nexttile;
plot(time, noise_real, 'Color', '#FF5757', 'LineWidth', 1.5);
xlim([0, n_symbols]); xlabel('símbolos');
title(sprintf('Ruido AWGN (SNR = %d dB)', SNR_dB));
grid on;

% Señal modulada con ruido
nexttile;
plot(time, noisy_signal, 'Color', '#23e623', 'LineWidth', 1.5);
xlim([0, n_symbols]); xlabel('símbolos');
title('Señal Modulada 32QAM con Ruido');
grid on;

%% Métricas de calidad
fprintf('=== MÉTRICAS DE CALIDAD DE LA SEÑAL ===\n');
fprintf('SNR configurada: %.1f dB\n', SNR_dB);
fprintf('Potencia de la señal: %.4f\n', signal_power);
fprintf('Potencia del ruido: %.4f\n', mean(abs(noise_real).^2));
fprintf('SNR medida: %.2f dB\n\n', 10*log10(signal_power/mean(abs(noise_real).^2)));

% Calcular Error Vector Magnitude (EVM)
evm = sqrt(mean(abs(noisy_symbols - qam_symbols).^2)) / sqrt(mean(abs(qam_symbols).^2));
fprintf('EVM (Error Vector Magnitude): %.2f%%\n', evm*100);
fprintf('Relación EVM a SNR: %.2f dB\n\n', -20*log10(evm));

%% Plot adicional: Zoom a una parte de la señal
figure('Position', [100 100 1000 400]);
tiledlayout(1, 2);

% Zoom señal sin ruido
nexttile;
zoom_start = floor(0.2 * length(time));
zoom_end = floor(0.25 * length(time));
plot(time(zoom_start:zoom_end), modulated_signal(zoom_start:zoom_end), 'Color', color, 'LineWidth', 2);
xlabel('Tiempo (s)');
title('Señal Modulada (sin ruido) - Zoom');
grid on;

% Zoom señal con ruido
nexttile;
plot(time(zoom_start:zoom_end), noisy_signal(zoom_start:zoom_end), 'Color', '#23e623', 'LineWidth', 2);
hold on;
plot(time(zoom_start:zoom_end), modulated_signal(zoom_start:zoom_end), 'Color', color, 'LineWidth', 1, 'LineStyle', '--');
hold off;
xlabel('Tiempo (s)');
title(sprintf('Señal con Ruido (SNR = %d dB) - Zoom', SNR_dB));
legend('Con ruido', 'Original', 'Location', 'best');
grid on;