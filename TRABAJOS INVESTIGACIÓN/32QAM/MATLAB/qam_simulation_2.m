clc; clear; close all;

%% Parámetros de la simulación
order = 32;
symorder = 'gray';
M = order;  % Número de símbolos
k = log2(M); % Bits por símbolo

% Parámetros de ruido (rango para análisis BER)
SNR_dB_values = 0:2:20;  % Valores de SNR a evaluar
num_SNR = length(SNR_dB_values);
cl
% Inicializar vectores para BER
BER_simulated = zeros(1, num_SNR);
BER_theoretical = zeros(1, num_SNR);

% Número de bits por simulación (debe ser múltiplo de k)
bits_per_simulation = 10000;
num_symbols_per_sim = bits_per_simulation / k;

%% Símbolos 32QAM ideales para referencia
symbols = 0 : M - 1;
bits_ref = reshape(de2bi(symbols).', [], 1);
qam_symbols_ref = qammod(bits_ref, M, symorder, "InputType", "bit", "UnitAveragePower", true);

%% Bucle principal para diferentes valores de SNR
for snr_idx = 1:num_SNR
    SNR_dB = SNR_dB_values(snr_idx);
    
    %% 1. GENERACIÓN DE BITS ALEATORIOS
    data_bits = randi([0 1], bits_per_simulation, 1);
    
    %% 2. MODULACIÓN 32QAM
    % Reorganizar bits en grupos de k (5 para 32QAM)
    reshaped_bits = reshape(data_bits, k, num_symbols_per_sim).';
    
    % Modulación
    tx_symbols = qammod(reshaped_bits, M, symorder, "InputType", "bit", "UnitAveragePower", true);
    
    %% 3. AGREGAR RUIDO AWGN
    % Calcular potencia de señal
    signal_power = mean(abs(tx_symbols).^2);
    
    % Convertir SNR de dB a lineal
    SNR_linear = 10^(SNR_dB/10);
    
    % Calcular potencia de ruido
    noise_power = signal_power / SNR_linear;
    
    % Generar ruido complejo
    noise = sqrt(noise_power/2) .* (randn(size(tx_symbols)) + 1i*randn(size(tx_symbols)));
    
    % Señal recibida con ruido
    rx_symbols = tx_symbols + noise;
    
    %% 4. DEMODULACIÓN 32QAM
    % Demodulación a bits
    rx_bits = qamdemod(rx_symbols, M, symorder, "OutputType", "bit", "UnitAveragePower", true);
    
    % Reorganizar bits a vector
    rx_bits = reshape(rx_bits.', [], 1);
    
    %% 5. CÁLCULO DE BER (Bit Error Rate)
    % Contar errores
    bit_errors = sum(data_bits ~= rx_bits);
    BER_simulated(snr_idx) = bit_errors / bits_per_simulation;
    
    %% 6. BER TEÓRICO PARA 32QAM
    % Para QAM rectangular, la probabilidad de error aproximada es:
    % P_b ≈ (4/k)*(1-1/sqrt(M))*Q(sqrt(3*k*SNR/(M-1)))
    
    % Función Q (complementary error function)
    Q = @(x) 0.5 * erfc(x/sqrt(2));
    
    % SNR por bit
    EbN0_linear = SNR_linear / k;
    
    % Probabilidad de error teórica (aproximación para QAM rectangular)
    % Esta es una aproximación común para QAM
    BER_theoretical(snr_idx) = (4/k) * (1 - 1/sqrt(M)) * Q(sqrt(3*k*EbN0_linear/(M-1)));
    
    %% 7. VISUALIZACIÓN PARA EL PRIMER Y ÚLTIMO VALOR DE SNR
    if snr_idx == 1 || snr_idx == num_SNR
        % Crear figura para constelación
        figure('Position', [100 100 1200 500]);
        
        % Constelación transmitida
        subplot(1, 2, 1);
        scatter(real(tx_symbols(1:min(1000, length(tx_symbols)))), ...
                imag(tx_symbols(1:min(1000, length(tx_symbols)))), ...
                30, 'b', 'filled');
        hold on;
        scatter(real(qam_symbols_ref), imag(qam_symbols_ref), 100, 'r', 'x', 'LineWidth', 2);
        hold off;
        xlim([-1.5 1.5]); ylim([-1.5 1.5]);
        grid on; axis square;
        title(sprintf('Constelación Transmitida (SNR = %d dB)', SNR_dB));
        xlabel('I Component'); ylabel('Q Component');
        legend('Símbolos TX', 'Posiciones Ideales', 'Location', 'best');
        
        % Constelación recibida
        subplot(1, 2, 2);
        scatter(real(rx_symbols(1:min(1000, length(rx_symbols)))), ...
                imag(rx_symbols(1:min(1000, length(rx_symbols)))), ...
                30, 'g', 'filled');
        hold on;
        scatter(real(qam_symbols_ref), imag(qam_symbols_ref), 100, 'r', 'x', 'LineWidth', 2);
        hold off;
        xlim([-1.5 1.5]); ylim([-1.5 1.5]);
        grid on; axis square;
        title(sprintf('Constelación Recibida (SNR = %d dB)', SNR_dB));
        xlabel('I Component'); ylabel('Q Component');
        legend('Símbolos RX', 'Posiciones Ideales', 'Location', 'best');
        
        % Calcular y mostrar métricas
        SER = sum(any(reshape(data_bits ~= rx_bits, k, []))) / num_symbols_per_sim;
        fprintf('SNR = %d dB:\n', SNR_dB);
        fprintf('  BER simulado: %.2e\n', BER_simulated(snr_idx));
        fprintf('  BER teórico:  %.2e\n', BER_theoretical(snr_idx));
        fprintf('  SER (Symbol Error Rate): %.2e\n', SER);
        fprintf('  Número de errores de bit: %d/%d\n\n', bit_errors, bits_per_simulation);
    end
end

%% 8. GRÁFICA COMPLETA DEL BER vs SNR
figure('Position', [100 100 900 600]);

% Gráfica semilogarítmica
semilogy(SNR_dB_values, BER_simulated, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', 'b', 'DisplayName', 'BER Simulado');
hold on;
semilogy(SNR_dB_values, BER_theoretical, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', 'r', 'DisplayName', 'BER Teórico (Aprox.)');

% Línea de referencia para BER = 1e-5
plot([min(SNR_dB_values) max(SNR_dB_values)], [1e-5 1e-5], 'k:', 'LineWidth', 1, ...
     'DisplayName', 'BER = 10^{-5}');

grid on;
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('Curva BER vs SNR para Modulación 32QAM', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
ylim([1e-7 1]);
xlim([min(SNR_dB_values) max(SNR_dB_values)]);

% Añadir anotaciones
text(mean(SNR_dB_values), 1e-2, '32QAM Modulation', ...
     'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(mean(SNR_dB_values), 5e-3, sprintf('M = %d, k = %d bits/símbolo', M, k), ...
     'FontSize', 10, 'HorizontalAlignment', 'center');

%% 9. DEMOSTRACIÓN DETALLADA PARA UN VALOR DE SNR ESPECÍFICO
SNR_demo = 15;  % SNR para demostración detallada
fprintf('\n=== DEMOSTRACIÓN DETALLADA PARA SNR = %d dB ===\n', SNR_demo);

% Generar nueva secuencia de bits para demostración
demo_bits = randi([0 1], 50, 1);  % Solo 50 bits para visualización clara
demo_symbols_num = length(demo_bits) / k;

% Modulación
demo_bits_reshaped = reshape(demo_bits, k, demo_symbols_num).';
demo_tx_symbols = qammod(demo_bits_reshaped, M, symorder, "InputType", "bit", "UnitAveragePower", true);

% Añadir ruido
demo_signal_power = mean(abs(demo_tx_symbols).^2);
demo_SNR_linear = 10^(SNR_demo/10);
demo_noise_power = demo_signal_power / demo_SNR_linear;
demo_noise = sqrt(demo_noise_power/2) .* (randn(size(demo_tx_symbols)) + 1i*randn(size(demo_tx_symbols)));
demo_rx_symbols = demo_tx_symbols + demo_noise;

% Demodulación
demo_rx_bits = qamdemod(demo_rx_symbols, M, symorder, "OutputType", "bit", "UnitAveragePower", true);
demo_rx_bits = reshape(demo_rx_bits.', [], 1);

% Calcular errores
demo_errors = demo_bits ~= demo_rx_bits;
demo_BER = sum(demo_errors) / length(demo_bits);

% Figura de demostración
figure('Position', [100 100 1400 800]);

% Subplot 1: Bits transmitidos vs recibidos
subplot(2, 2, 1);
stem(1:length(demo_bits), demo_bits, 'b', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'b');
hold on;
stem(find(demo_errors), demo_bits(demo_errors), 'r', 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 10);
stem(1:length(demo_rx_bits), demo_rx_bits, 'g--', 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 4);
hold off;
xlabel('Índice de Bit', 'FontSize', 11);
ylabel('Valor del Bit', 'FontSize', 11);
title(sprintf('Bits Transmitidos vs Recibidos (SNR = %d dB)', SNR_demo), 'FontSize', 12);
legend('Bits TX', 'Errores', 'Bits RX', 'Location', 'best');
grid on;
xlim([0 length(demo_bits)+1]);
ylim([-0.2 1.2]);

% Subplot 2: Constelación de demostración
subplot(2, 2, 2);
scatter(real(demo_rx_symbols), imag(demo_rx_symbols), 80, 'g', 'filled', 'MarkerEdgeColor', 'k');
hold on;
scatter(real(demo_tx_symbols), imag(demo_tx_symbols), 120, 'b', 'x', 'LineWidth', 2);
scatter(real(qam_symbols_ref), imag(qam_symbols_ref), 200, 'r', '+', 'LineWidth', 2);
hold off;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
grid on; axis square;
xlabel('Componente I', 'FontSize', 11);
ylabel('Componente Q', 'FontSize', 11);
title('Constelación: TX vs RX vs Ideal', 'FontSize', 12);
legend('RX (con ruido)', 'TX (original)', 'Posiciones Ideales', 'Location', 'best');

% Subplot 3: Distribución de errores por símbolo
subplot(2, 2, 3);
symbol_errors = any(reshape(demo_errors, k, []));
symbol_error_positions = find(symbol_errors);
bar(1:demo_symbols_num, double(symbol_errors), 'FaceColor', 'r', 'EdgeColor', 'k');
xlabel('Índice de Símbolo', 'FontSize', 11);
ylabel('Error (0/1)', 'FontSize', 11);
title(sprintf('Errores por Símbolo (BER = %.3f)', demo_BER), 'FontSize', 12);
grid on;
xlim([0 demo_symbols_num+1]);
ylim([0 1.2]);

% Subplot 4: Comparación BER teórico vs simulado
subplot(2, 2, 4);
bar_values = [BER_theoretical(SNR_dB_values == SNR_demo), demo_BER];
bar(1:2, bar_values, 'FaceColor', [0.5 0.5 0.8]);
set(gca, 'XTickLabel', {'Teórico', 'Simulado'});
ylabel('BER', 'FontSize', 11);
title('Comparación BER Teórico vs Simulado', 'FontSize', 12);
text(1:1, bar_values, num2str(bar_values', '%.2e'), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
grid on;
ylim([0 max(bar_values)*1.5]);

%% 10. TABLA RESUMEN DE RESULTADOS
fprintf('\n=== TABLA RESUMEN DE RESULTADOS BER ===\n');
fprintf('SNR (dB)\tBER Simulado\tBER Teórico\tDiferencia\n');
fprintf('----------------------------------------------------\n');
for i = 1:num_SNR
    diff = abs(BER_simulated(i) - BER_theoretical(i));
    fprintf('%d\t\t%.2e\t\t%.2e\t\t%.2e\n', ...
            SNR_dB_values(i), BER_simulated(i), BER_theoretical(i), diff);
end

%% 11. ANÁLISIS ADICIONAL: EVM vs SNR
fprintf('\n=== ANÁLISIS DE EVM (Error Vector Magnitude) ===\n');

% Recalcular EVM para cada SNR
EVM_values = zeros(1, num_SNR);

for snr_idx = 1:num_SNR
    SNR_dB = SNR_dB_values(snr_idx);
    
    % Generar nueva señal
    test_bits = randi([0 1], 10000, 1);
    test_symbols = qammod(reshape(test_bits, k, []).', M, symorder, "InputType", "bit", "UnitAveragePower", true);
    
    % Añadir ruido
    test_power = mean(abs(test_symbols).^2);
    test_SNR_linear = 10^(SNR_dB/10);
    test_noise = sqrt(test_power/(2*test_SNR_linear)) .* (randn(size(test_symbols)) + 1i*randn(size(test_symbols)));
    test_rx_symbols = test_symbols + test_noise;
    
    % Calcular EVM
    EVM_values(snr_idx) = sqrt(mean(abs(test_rx_symbols - test_symbols).^2)) / sqrt(mean(abs(test_symbols).^2));
end

% Gráfica EVM vs SNR
figure('Position', [100 100 900 400]);

subplot(1, 2, 1);
plot(SNR_dB_values, 100*EVM_values, 'm-^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
xlabel('SNR (dB)', 'FontSize', 11);
ylabel('EVM (%)', 'FontSize', 11);
title('EVM vs SNR para 32QAM', 'FontSize', 12);
grid on;

subplot(1, 2, 2);
plot(SNR_dB_values, 20*log10(EVM_values), 'c-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'c');
xlabel('SNR (dB)', 'FontSize', 11);
ylabel('EVM (dB)', 'FontSize', 11);
title('EVM en dB vs SNR', 'FontSize', 12);
grid on;

% Mostrar relación EVM-SNR
fprintf('\nRelación aproximada: EVM (dB) ≈ -SNR (dB)\n');
fprintf('Para SNR alto, EVM (%%RMS) ≈ 100 * 10^{-SNR/20}\n');