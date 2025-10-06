clc
clear all
%signal conditioning and muxing for 3 input

[u1 fs1]=audioread('voz_1.au');
[u2 fs2]=audioread('voz_2.au');
[u3 fs3]=audioread('voz_3.au');

N=3;
ts=1/fs1;
Ts = ts; % Slot time interval
ns = Ts/ts;

for s=1:N    %TO find LENGTH
    eval(['zz'  '=u' num2str(s)])
    eval(['l' num2str(s)  '=length(zz)'])
end

% m1=max(l1,l2);
% m=max(m1,l3);
mm=0;

for s=1:N    % FOR MAX LENGTH
    if eval(['l' num2str(s)])>mm
        mm=eval(['l' num2str(s)]);
    end
    
end
m=mm;

for s=1:N    %TO MAKE EQUAL LENGTH

    if eval(['l' num2str(s) ' <m;'])
        eval(['ll' '=l' num2str(s) ])

        eval(['zz' '=u' num2str(s)  ])
        eval(['uu' num2str(s)  '=zz;'])   
        for i=ll+1:m
            eval(['uu' num2str(s)  '(i)=0;'])
        end
    else
        eval(['zz' '=u' num2str(s)  ])
        eval(['uu' num2str(s)  '=zz;']) 
    end
end
clc

disp('Press ENTER to hear the audio signals before conditioning');
pause
sound(uu1)
disp('Press ENTER to hear the audio signals before conditioning');
pause
sound(uu2)
disp('Press ENTER to hear the audio signals before conditioning');
pause
sound(uu3)


for b=1:N
    eval(['u(b,:)=uu' num2str(b) ';' ])
end 


uuu=TDM_mux2(u,ns);

disp('Hit [Enter}]to listen to the muxed sound');
pause
sound(uuu)

y=demux_TDM2(uuu,N,ns);
