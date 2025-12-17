%% Iniciando variables
close all;clear all; clc;

% Acondicionamiento de la señales para la multiplexacion para N entradas

N=3; % numero de canales (fuentes de informacion)

for a=1:N
    if(N<10)
        eval(['[u' num2str(a) ' fs' num2str(a) ']=audioread(''voz_' num2str(a) '.au'');' ])
    else
        eval(['[u' num2str(a) ' fs' num2str(a) ']=audioread(''voz_' num2str(a) '.au'');' ])
    end
end % loop end

% periodo de muestreo
ts=1/fs1; 
% Duracion del slot
Ts = ts; % 

ns = Ts/ts;

% Escuchar las señales de voz
disp('Presiona enter para continuar');
pause
disp('Presiona entre para reproducir el primer archivo de audio');
sound(u1)
pause
disp('Presiona entre para reproducir el segundo archivo de audio');
sound(u2)
pause
disp('Presiona entre para reproducir el tercer archivo de audio');
sound(u3)

%% Ajustes
% Determinar la longitud de las muestras
for s=1:N    
    eval(['zz'  '=u' num2str(s) ';'])
    eval(['l' num2str(s)  '=length(zz);'])
end

mm=0;

% Determinar la máxima longitud
for s=1:N    
    if eval(['l' num2str(s)])>mm
        mm=eval(['l' num2str(s) ';']);
    end    
end
m=mm;

% Estimacion de tiempo
tm=ts*(m-1);
tx=0:ts:tm;

% Igualando la longitud de los vectores
for s=1:N

    if eval(['l' num2str(s) ' <m'])
        eval(['ll' '=l' num2str(s) ';' ]);
        eval(['zz' '=u' num2str(s) ';'  ])
        eval(['uu' num2str(s)  '=zz;'])  ;
        eval(['uu' num2str(s)  '(ll+1:m)=0;']);
    else
        eval(['zz' '=u' num2str(s) ';']);
        eval(['uu' num2str(s)  '=zz;']) ;
    end
end

%% Multiplexando (union en una matriz)
for b=1:N
    eval(['u(' num2str(b) ',:)=uu' num2str(b) ';' ])
end % loop end


% Graficando las senales
figure(1)

for i=1:N
    subplot(N,1,i)
    stairs(tx,u(i,:));
    xlabel('Tiempo [s]');
    ylabel('Amplitud');  
    if i==1
        title('Senal de voz')
    end
end

%% Procesamiento
% Multiplexando la senal

um=tdmm(u,ns);
% Agregando AWGN (ruido gaussiano aditivo aleatorio)
umn=awgn(um,25);

% Plotting: TDM signal with and w/o noise
figure(2)
subplot(2,1,1)

stairs(tx,um);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Senal TDM')

subplot(2,1,2)
stairs(tx,umn);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Senal + ruido')


disp('Presiona enter para escuchar la señal multiplexada: ');
pause
sound(um)

% Procesando (operacion inversa) - Demultiplexando
% um -> [Demux] -> ud
ud=demux(um,ns,N);

% Dividiendo las senales demultiplexadas
for b=1:N
    eval(['ux' num2str(b) '=ud(' num2str(b) ',:);' ]);
end

% Reproduciendo las senales demultiplexadas 
display('Presiona enter para escuchar las senales demultiplexadas')
pause
for b=1:N
    
    eval(['sound(ux' num2str(b) ');' ])
    display('Presionar enter para reproducir la senal demultiplexada')
    pause
    
end % loop end

% Graficando las senales demultiplexadas
figure(3)
for a=1:N
    subplot(N,1,a)
    stairs(tx,ud(a,:));
    xlabel('Tiempo [s]');
    ylabel('Amplitud'); 
    if a==1
        title('Senal Demultiplexada')
    end
end

