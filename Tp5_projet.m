close all
clear all
clc

T = 10;
fs = 1000;
t = 0:1/fs:T;
N = length(t);

f1 = 50;
f2 = 100;

A1 = 2;
A2 = 1;

% Signal combiné
s1 = A1*sin(2*pi*f1*t);
s2 = A2*sin(2*pi*f2*t);
signal = s1 + s2;
subplot(2,1,1)
plot(t,signal)
fshift = (-N/2:N/2-1)*(fs/N);
f = (0:N-1)*(fs/N);
y = fft(signal);
subplot(2,1,2)
plot(fshift, fftshift(abs(y)));

sound(s1);
%%

%creation des filtres

%fonction Butter

%filtre pass_bas
fc_bas=f2/(fs/2);%la fréquence de coupure pour faire passer la fréquence du signal 2 (100) et éliminer la fréquence (50)
[b_low,a_low] = butter(1,fc_bas,'low');


%filtre_pass_haut
fc_haut=f1/(fs/2);%la fréquence de coupure pour faire passer la fréquence du signal 1 (50) et éliminer la fréquence (100)
[b_high,a_high] = butter(4,fc_haut,'high');


%tracage des signaux 
figure(1)
freqz(b_low,a_low);
title("la réponse fréquentielle du filtre numérique pass_bas ");
grid on
figure(2)
freqz(b_high,a_high);
title("la réponse fréquentielle du filtre numérique pass_haut ");
grid on

%%
%{
%fonction Cheby

%filtre pass_bas
%fc_bas=f2/(fs/2);%la fréquence de coupure pour faire passer la fréquence du signal 2 (100) et éliminer la fréquence (50)
[b1_low,a1_low] = cheby1(4,fc_bas,'low');


%filtre_pass_haut
%fc_haut=f1/(fs/2);%la fréquence de coupure pour faire passer la fréquence du signal 1 (50) et éliminer la fréquence (100)
[b1_high,a1_high] = cheby1(4,fc_haut,'high');


%tracage des signaux 
figure(3)
freqz(b1_low,a1_low);
title("la réponse fréquentielle du filtre numérique pass_bas ");
grid on
figure(4)
freqz(b1_high,a1_high);
title("la réponse fréquentielle du filtre numérique pass_haut ");
grid on

%}

%% 
%les signaux de sorties :
s1 = A1*sin(2*pi*f1*t);
s2 = A2*sin(2*pi*f2*t);
signal = s1 + s2;
fc_bas=f2/(fs/2);%la fréquence de coupure pour faire passer la fréquence du signal 2 (100) et éliminer la fréquence (50)
[b_low,a_low] = butter(1,fc_bas,'low');
y_low = filter(b_low, a_low, signal);
sound(y_low);
%% filtrage de sons 

[son1,Fe1]=audioread("affreux.wav");
[son2,Fe2]=audioread("akesuper.wav");

N1=length(son1);
N2=length(son2);

Te1 = 1/Fe1;
t1 = 0:Te1:(N1-1)*Te1;

Te2 = 1/Fe2;
t2 = 0:Te2:(N2-1)*Te2;


subplot(221)
plot(son1);
title('representation temporelle du signal audio1 ')
xlabel('t')
ylabel('Amplitude')
grid on

subplot(223)
plot(son2);
title('representation temporelle du signal audio2 ')
xlabel('t')
ylabel('Amplitude')
grid on

sound(son1)

%sound(son2)

%%Fe1_2=Fe1*1.15;%la frequence d'echantillonage normale
%%Fe2_2=Fe2*1.15;

%Analyse frequentielle du son1

fshift = (-N1/2:N1/2-1)*(Fe1/N1); % le pas de discrétisation : Fe/N
transfF=fft(son1); % transformée de fourirer

subplot(222)
plot(fshift,fftshift(abs(transfF)));
title('representation fréquentielle du signal audio1 ')
xlabel('f')
ylabel('Amplitude')
grid on

%Analyse frequentielle du son2

fshift2 = (-N2/2:N2/2-1)*(Fe2/N2); % le pas de discrétisation : Fe/N
transfF2=fft(son2); % transformée de fourirer

subplot(224)
plot(fshift2,fftshift(abs(transfF2)));
title('representation fréquentielle du signal audio2 ')
xlabel('f')
ylabel('Amplitude')
grid on

%%

subplot(221)
snoise1 = son1+0.05*randn(size(son1));
plot(t1,snoise1)
xlabel('t');
ylabel('snoise1');
title('Représentation temporelle du signal bruité');
grid on


subplot(222)
ysnoise1=fft(snoise1);
plot(fshift,fftshift(abs(ysnoise1)*2/N1));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');
grid on

subplot(223)
snoise2 = son2+0.05*randn(size(son2));
plot(t2,snoise2)
xlabel('t');
ylabel('snoise2');
title('Représentation temporelle du signal bruité');
grid on

subplot(224)
ysnoise2=fft(snoise2);
plot(fshift2,fftshift(abs(ysnoise2)*2/N2));
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal bruité');
grid on

sound(snoise1)

%% filtrage numerique

[a_noise1,b_noise1] = butter(1,0.6,"low")


figure(1)
freqz(a_noise1,b_noise1);
title("la réponse fréquentielle du filtre numérique pass_bas ");
grid on

son1_filter = filter(a_noise1,b_noise1,snoise1);
sound(son1_filter)

%%[a_noise2,b_noise2] = butter(1,0.4,"high")
%%son2_filter = filter(a_noise1,b_noise2,son1_filter);
%sound(son2_filter)

% son1_filter_FFT = fft(son1_filter);
% plot(fshift,fftshift(abs(son1_filter)));

% Enregistrement du signal filtré
wavwrite(son1_filter,Fe1,'bruit_filtre.wav');