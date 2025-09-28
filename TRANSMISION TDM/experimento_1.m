% TUTORIAL / EJEMPLO MATLAB: Generación de señales, PCM, TDM y recuperación
% Dirigido a estudiantes de telecomunicaciones (nivel intermedio).
% -------------------------------------------------------------
clear; close all; clc;

%% Parámetros de señal continua (simulación "alta resolución")
Fs_cont = 20000;          % frecuencia de muestreo "continua" para simular analogico (Hz)
T_cont = 0.1;             % duración de la señal (s)
t_cont = 0:1/Fs_cont:T_cont-1/Fs_cont;   % vector tiempo continuo

% Parámetros de cada señal (elige frecuencias que representen una banda base)
f_seno = 120;             % frecuencia seno (Hz)
f_tri  = 60;              % frecuencia triangular (Hz)
A_seno = 2;             % amplitud seno
A_tri  = 2;             % amplitud triangular

%% 1) Generación de la señal senoidal y triangular
x_seno = A_seno * sin(2*pi*f_seno*t_cont);

% Triangular sin toolboxes: señal triangular (amplitud entre -1 y 1)
% fórmula: tri(t) = 2 * abs( 2*(f*t - floor(f*t + 0.5)) ) - 1
x_tri = A_tri * (2*abs(2*(f_tri*t_cont - floor(f_tri*t_cont + 0.5))) - 1);

%% 2) Decidir frecuencia de muestreo para muestreador (respetando Nyquist)
% Teorema de Nyquist: fs >= 2 * f_max (f_max = la componente más alta)
f_max = max([f_seno, f_tri]);
% Elegimos sobremuestreo razonable (factor 10 para reducir aliasing en ejemplo)
factor_sobre = 10;
fs = max(2*f_max, factor_sobre * f_max);
fs = ceil(fs);            % frecuencia de muestreo usada por el muestreador (Hz)
Ts = 1/fs;

% Tiempo muestreado (vector de instantes de muestreo)
t_sample = 0:Ts:T_cont-Ts;
Nsamples = length(t_sample);

%% 3) Muestreo (emulación del muestreador)
% Muestreamos las señales continuas en los instantes t_sample (muestreo ideal)
x_seno_sampled = A_seno * sin(2*pi*f_seno*t_sample);
x_tri_sampled  = A_tri  * (2*abs(2*(f_tri*t_sample - floor(f_tri*t_sample + 0.5))) - 1);

% Para visualizaciones, también guardamos la versión muestreada como impulsiva
% (pero las señales sampleadas arriba sirven para PCM).

%% 4) Conversión a PCM (cuantización uniforme y codificación binaria)
Nbits = 8;                        % número de bits por muestra (puedes cambiarlo)
Nlevels = 2^Nbits;

% Normalizamos cada señal al rango [-1,1] antes de cuantizar (evitar saturaciones).
% Si ya conoces el rango real usalo en lugar de normalizar (aquí usamos amplitud maxima)
% Normalización por el valor máximo absoluto de cada señal muestreada:
max_seno = max(abs(x_seno_sampled));
max_tri  = max(abs(x_tri_sampled));
x_seno_norm = x_seno_sampled / max_seno;
x_tri_norm  = x_tri_sampled  / max_tri;

% Cuantizador uniforme: niveles de - (2^(N-1))..(2^(N-1)-1) (representación signed)
Qmax = 2^(Nbits-1) - 1;
% cuantización (enteros)
q_seno = round( x_seno_norm * Qmax );
q_tri  = round( x_tri_norm  * Qmax );

% Asegurar que no se desborde
q_seno(q_seno > Qmax) = Qmax;
q_seno(q_seno < -Qmax) = -Qmax;
q_tri(q_tri > Qmax) = Qmax;
q_tri(q_tri < -Qmax) = -Qmax;

% Reconstrucción cuantizada (nivel cuantizado a amplitud real)
x_seno_cuant = (q_seno / Qmax) * max_seno;
x_tri_cuant  = (q_tri  / Qmax) * max_tri;

% Codificación binaria (PCM): convertimos cada muestra (signed) a palabra Nbits
% Usamos representación en complemento a dos para signos.
% Para convertir en complemento a dos a cadena binaria trabajamos con uint16 offset.
offset = 2^(Nbits-1); % para pasar a un rango 0..2^N-1 (si usamos complemento a dos)
q_seno_unsigned = uint16( q_seno + offset ); % ahora 0..2^N-1
q_tri_unsigned  = uint16( q_tri  + offset );

pcm_seno_bin = dec2bin(q_seno_unsigned, Nbits); % matriz char: Nsamples x Nbits
pcm_tri_bin  = dec2bin(q_tri_unsigned, Nbits);

% Crear un bitstream (secuencial) por canal (concatenando palabras)
% (cada fila es una muestra; para bitstream por canal concatenamos filas)
pcm_seno_bitstream = reshape(pcm_seno_bin', 1, []); % char vector con '0'/'1'
pcm_tri_bitstream  = reshape(pcm_tri_bin',  1, []);

%% 5) Emulación TDM
% ---- 5a) TDM en dominio de muestras (intercalado de muestras)
% Estructura de trama simple: por cada instante n, enviamos muestra Seno(n) luego Tri(n)
% Esto es TDM a nivel de muestras (no bits).
tdm_muestras = zeros(1, 2*Nsamples);
tdm_muestras(1:2:end) = x_seno_sampled;   % posiciones impares
tdm_muestras(2:2:end) = x_tri_sampled;    % posiciones pares

% Tiempo asociado a la señal TDM (cada muestra ocupa Ts/2 en este esquema)
Ts_tdm = Ts / 2;
t_tdm = 0:Ts_tdm:(length(tdm_muestras)-1)*Ts_tdm;

% ---- 5b) TDM de bits: trama concatenada intercalando palabras binarias
% Construimos bitstream TDM alternando palabras completas de cada canal
% Ejemplo: [seno_word1, tri_word1, seno_word2, tri_word2, ...]
% Preasignar tamaño: cada muestra aporta Nbits
tdm_bits_len = 2 * Nsamples * Nbits;
tdm_bitstream = char(zeros(1, tdm_bits_len));
idx = 1;
for n = 1:Nsamples
    % palabra seno
    tdm_bitstream(idx:idx+Nbits-1) = pcm_seno_bin(n,:);
    idx = idx + Nbits;
    % palabra tri
    tdm_bitstream(idx:idx+Nbits-1) = pcm_tri_bin(n,:);
    idx = idx + Nbits;
end

%% 6) Recuperación (demultiplexado) y reconstrucción de señales
% Demultiplexado en dominio de muestras
rec_seno_samples = tdm_muestras(1:2:end); % tomar posiciones impares
rec_tri_samples  = tdm_muestras(2:2:end); % posiciones pares

% Reconstrucción en tiempo continuo (interpolación). Usamos interp1 con 'sinc-like' (pchip o spline)
% Vector de tiempos muestreados originales:
t_rec = t_sample; % Ts muestreo original
% Queremos reconstruir al tiempo de alta resolución t_cont
% Usaremos interpolación 'spline' para aproximar una reconstrucción (no es ideal, pero instructiva)
x_seno_rec = interp1(t_rec, rec_seno_samples, t_cont, 'spline');
x_tri_rec  = interp1(t_rec, rec_tri_samples,  t_cont, 'spline');

% Alternativa: reconstrucción por hold (ZOH) para mostrar diferencia:
x_seno_zoh = zeros(size(t_cont));
x_tri_zoh = zeros(size(t_cont));
for k=1:Nsamples
    idxs = (t_cont >= t_rec(k)) & (t_cont < t_rec(min(k+1,Nsamples)));
    x_seno_zoh(idxs) = rec_seno_samples(k);
    x_tri_zoh(idxs)  = rec_tri_samples(k);
end

%% 7) Gráficas (todas en una figura con subplots)
figure('Name','TDM y PCM - Ejemplo','Units','normalized','Position',[0.05 0.05 0.9 0.85]);

% 1: Señal original (alta resolución)
subplot(3,2,1);
plot(t_cont, x_seno, 'LineWidth', 1);
title('Señal original: Seno');
xlabel('Tiempo (s)'); ylabel('Amplitud');
grid on;

subplot(3,2,2);
plot(t_cont, x_tri, 'LineWidth', 1);
title('Señal original: Triangular');
xlabel('Tiempo (s)'); ylabel('Amplitud');
grid on;

% 2: Señales muestreadas y cuantizadas (muestras marcadas)
subplot(3,2,3);
plot(t_cont, x_seno, 'LineWidth', 0.8); hold on;
stem(t_sample, x_seno_sampled, 'filled'); % muestras
plot(t_sample, x_seno_cuant, 'o','MarkerSize',4); % muestras cuantizadas
title(sprintf('Seno: Muestreada y cuantizada (Nbits=%d)',Nbits));
xlabel('Tiempo (s)'); ylabel('Amplitud'); legend('Original','Muestras','Cuantizadas');
grid on; hold off;

subplot(3,2,4);
plot(t_cont, x_tri, 'LineWidth', 0.8); hold on;
stem(t_sample, x_tri_sampled, 'filled');
plot(t_sample, x_tri_cuant, 'o','MarkerSize',4);
title(sprintf('Tri: Muestreada y cuantizada (Nbits=%d)',Nbits));
xlabel('Tiempo (s)'); ylabel('Amplitud'); legend('Original','Muestras','Cuantizadas');
grid on; hold off;

% 3: TDM (muestras intercaladas) — mostramos fragmento para claridad
subplot(3,2,5);
stairs(t_tdm, tdm_muestras, 'LineWidth', 1);
xlim([0, 5*Ts_tdm]); % mostrar unas pocas tramas
title('TDM (muestras intercaladas) - fragmento');
xlabel('Tiempo (s)'); ylabel('Amplitud'); grid on;

% 4: Señales recuperadas (comparación con original)
subplot(3,2,6);
plot(t_cont, x_seno, 'LineWidth', 0.8); hold on;
plot(t_cont, x_seno_rec, '--','LineWidth', 1.2);
plot(t_cont, x_seno_zoh, ':','LineWidth', 1);
plot(t_cont, x_tri, 'LineWidth', 0.6);
plot(t_cont, x_tri_rec, '--','LineWidth', 1.2);
legend('Seno orig','Seno rec spline','Seno ZOH','Tri orig','Tri rec spline');
title('Señales originales y recuperadas (spline y ZOH)');
xlabel('Tiempo (s)'); ylabel('Amplitud'); grid on; hold off;

% Ajustes finales: mejorar espaciado
sgtitle('Flujo: Generación -> Muestreo -> PCM -> TDM -> Recuperación');

%% 8) Mostrar métricas (SNR de cuantización aproximada)
% error cuantización por canal
err_seno = x_seno_sampled - x_seno_cuant;
err_tri  = x_tri_sampled  - x_tri_cuant;
SNR_seno_db = 10*log10( sum(x_seno_sampled.^2) / sum(err_seno.^2) );
SNR_tri_db  = 10*log10( sum(x_tri_sampled.^2)  / sum(err_tri.^2) );

fprintf('Frecuencia de muestreo usada fs = %d Hz (Ts = %.5f s)\n', fs, Ts);
fprintf('Muestras por canal = %d, Bits por muestra = %d\n', Nsamples, Nbits);
fprintf('SNR aproximada cuantizacion: Seno = %.2f dB, Tri = %.2f dB\n', SNR_seno_db, SNR_tri_db);

%% FIN DEL SCRIPT

% -------------------------------------------------------------------------
% En caso de requerir más especificaciones para ejecutar lo solicitado:
% - Indicar valores concretos de frecuencia y amplitud si no desea los de ejemplo.
% - Indicar Nbits deseados (por defecto Nbits=8).
% - Indicar si desea que la reconstrucción use un filtro ideal (sinc) en vez de 'spline'
%   (la implementación de sinc requiere ventana y convolución; con gusto la añado).
% - Indicar si quiere que la TDM sea por trama de tiempo fija con preámbulo, sincronismo, CRC, etc.
% - Indicar si desea que la salida TDM se exporte a un archivo binario (.bin) o a WAV (para audio).
% -------------------------------------------------------------------------
