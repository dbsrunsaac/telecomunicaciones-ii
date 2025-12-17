clc;
state = 1;
% Solo realiza la grabación cuando el estado es 0
if state == 1 
    clear;
    % Objeto de grabacion
    recObj = audiorecorder;
    
    % Iniciar una grabacion (solo es necesario ejecutarla una vez por audio)
    disp("Iniciando con la grabacion");
    recordblocking(recObj, 3);
end

% Reproducir el audio grabado
play(recObj)
y = getaudiodata(recObj);

save_voice = 1;

if save_voice == 1
    % guardar la señal de voz
    fs = 8000;
    audiowrite('voz_3.wav', y, fs);
end









