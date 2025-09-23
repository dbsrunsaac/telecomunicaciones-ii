% Código en MATLAB para generar señales, convertir a PCM, muestrear, multiplexar en TDM y recuperar señales.
% Este código está diseñado para estudiantes de mitad de carrera en telecomunicaciones, con explicaciones detalladas.
% Nivel intermedio: Utiliza funciones básicas de MATLAB para señales (sin toolboxes especiales como Signal Processing Toolbox,
% aunque se asume que está disponible para funciones como sawtooth y stairs).

% Paso 1: Generación de una señal triangular y senoidal.
% - Señal triangular: Utilizamos la función sawtooth(t, width) donde width=0.5 genera una onda triangular simétrica.
%   La frecuencia se define como f_tri, y el tiempo t se genera con linspace para un vector uniforme.
% - Señal senoidal: Utilizamos sin(2*pi*f_sin*t) para una onda sinusoidal estándar.
% Definimos parámetros iniciales:
f_tri = 1;      % Frecuencia de la señal triangular en Hz (baja para visualización clara).
f_sin = 2;      % Frecuencia de la señal senoidal en Hz (diferente para distinguir).
A = 1;          % Amplitud de ambas señales (normalizada para simplicidad).
duracion = 2;   % Duración total de la señal en segundos.
fs_cont = 1000; % Frecuencia de muestreo alta para simular señal continua (mucho mayor que Nyquist para plotting suave).
t = linspace(0, duracion, duracion * fs_cont);  % Vector de tiempo continuo usando linspace(num_inicio, num_fin, num_puntos).

% Generar señales:
senal_tri = A * sawtooth(2*pi*f_tri*t, 0.5);  % sawtooth genera onda diente de sierra; con 0.5 es triangular.
senal_sin = A * sin(2*pi*f_sin*t);            % sin() para sinusoidal.

% Plotear señales originales para verificación (opcional, pero útil para entender).
figure(1);
subplot(2,1,1);
plot(t, senal_tri);
title('Señal Triangular Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,1,2);
plot(t, senal_sin);
title('Señal Senoidal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Paso 2: Conversión de las señales a PCM.
% PCM (Pulse Code Modulation) involucra: muestreo (ver paso 3), cuantización y codificación.
% Aquí, preparamos para cuantización asumiendo muestreo posterior.
% Definimos parámetros para PCM:
n_bits = 8;     % Número de bits para cuantización (resolución 2^n niveles).
V_max = A;      % Voltaje máximo (amplitud máxima de la señal).
V_min = -A;     % Voltaje mínimo (para señales bipolares como estas).
niveles = 2^n_bits;  % Número de niveles de cuantización.

% La cuantización se hará después del muestreo en el paso 3, ya que PCM requiere muestreo primero.
% Función para cuantizar: Usaremos round() para cuantización uniforme.
delta = (V_max - V_min) / (niveles - 1);  % Paso de cuantización.

% Paso 3: Emular un muestreador con frecuencia de muestreo adecuada para Nyquist.
% Teorema de Nyquist: fs >= 2 * f_max para evitar aliasing.
% f_max = max(f_tri, f_sin) = 2 Hz, así que fs >= 4 Hz. Usamos fs = 10 Hz para margen.
fs = 10;  % Frecuencia de muestreo en Hz.
Ts = 1/fs;  % Período de muestreo.
num_muestras = floor(duracion / Ts) + 1;  % Número de muestras (floor para entero).
t_muestras = 0:Ts:duracion;  % Vector de tiempo muestreado (usando colon operator para pasos uniformes).

% Muestrear señales: Interpolamos las señales continuas en los puntos de muestreo.
% Usamos interp1 para interpolar (lineal por defecto, adecuado para señales suaves).
senal_tri_muest = interp1(t, senal_tri, t_muestras, 'linear');
senal_sin_muest = interp1(t, senal_sin, t_muestras, 'linear');

% Ahora, aplicar cuantización para completar PCM.
% Cuantizar: Redondear a niveles más cercanos.
senal_tri_quant = V_min + round((senal_tri_muest - V_min) / delta) * delta;
senal_sin_quant = V_min + round((senal_sin_muest - V_min) / delta) * delta;

% Codificación: Convertir a binario (representación como strings para simplicidad, no esencial para TDM).
% Para TDM, usaremos las muestras cuantizadas directamente (valores numéricos), ya que TDM multiplexa en tiempo.

% Paso 4: Obtener la señal TDM a partir del muestreador para cada señal.
% TDM (Time Division Multiplexing): Para 2 señales, intercalamos muestras: [muestra1_tri, muestra1_sin, muestra2_tri, muestra2_sin, ...]
% Asumimos multiplexación síncrona simple (mismo fs para ambas).
% Vector TDM: Usamos reshape y transpose para intercalar.
senal_tdm = zeros(1, 2 * num_muestras);  % Prealocar vector TDM (doble longitud).
senal_tdm(1:2:end) = senal_tri_quant;    % Posiciones impares: señal triangular (1-based indexing en MATLAB).
senal_tdm(2:2:end) = senal_sin_quant;    % Posiciones pares: señal senoidal.

% Tiempo para TDM: Cada "slot" es Ts/2, ya que multiplexamos 2 señales.
Ts_tdm = Ts / 2;  % Período efectivo por slot.
t_tdm = 0:Ts_tdm:(2*num_muestras-1)*Ts_tdm;  % Vector de tiempo para TDM.

% Plotear TDM para verificación.
figure(2);
stairs(t_tdm, senal_tdm);  % stairs() para plotear como escalones, simulando señal digital multiplexada.
title('Señal TDM Multiplexada');
xlabel('Tiempo (s)');
ylabel('Amplitud Cuantizada');
grid on;

% Paso 5: Recuperar cada señal transmitida por TDM y mostrar en una figura con subplot.
% Demultiplexar: Extraer impares para tri, pares para sin.
senal_tri_rec = senal_tdm(1:2:end);  % Recuperar triangular.
senal_sin_rec = senal_tdm(2:2:end);  % Recuperar senoidal.

% Para reconstruir, usamos hold (zero-order hold) simulando con stairs, o interpolamos de vuelta.
% Aquí, ploteamos las recuperadas en tiempo original de muestreo.

% Figura con subplots: Originales vs Recuperadas.
figure(3);
subplot(2,2,1);
plot(t, senal_tri);
title('Triangular Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,2,2);
stairs(t_muestras, senal_tri_rec);  % Usamos stairs para mostrar como muestras hold.
title('Triangular Recuperada de TDM');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,2,3);
plot(t, senal_sin);
title('Senoidal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(2,2,4);
stairs(t_muestras, senal_sin_rec);
title('Senoidal Recuperada de TDM');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Notas adicionales:
% - En un sistema real, TDM incluiría sincronización (e.g., frames), pero aquí es simplificado.
% - No se incluye ruido o transmisión; asume canal ideal.
% - Para reconstrucción perfecta, usar filtro paso bajo post-demux, pero omitido por simplicidad.

% Especificaciones adicionales requeridas:
% - Este código asume MATLAB base con funciones como sawtooth (de Signal Processing Toolbox). Si no está disponible, reemplaza sawtooth con una implementación manual: senal_tri = 2*A*(t*f_tri - floor(t*f_tri + 0.5)).
% - Ajusta fs, f_tri, f_sin según necesidades (e.g., para aliasing, reduce fs < 4 Hz).
% - Para más precisión en PCM, considera mu-law o A-law, pero aquí es lineal.
% - Ejecuta en MATLAB R2020 o superior para compatibilidad.