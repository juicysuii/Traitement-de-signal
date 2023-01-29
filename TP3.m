close all
clear all
clc

%Représentation temporelle

load('ecg.mat');
fe = 500; %fréquence
te = 1/fe;
N = length(ecg);
t = 0:te:(N-1)*te; %intervalle

plot(t,ecg);
grid on
xlabel('t')
ylabel('ECG')
title("Représentation temporelle de l'activation électrique du cœur")
xlim([0.5 1.5]);%zoom

%Représentation fréquentielle
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*fe/N;
y = fft(ecg);
plot(fshift,fftshift(abs(y)))
grid on
xlabel('f')
ylabel('Amplitude')
title("Représentation fréquentielle de l'activation électrique du cœur")

%Filtre passe_haut

filtre_ps_haut = ones(size(ecg));
fc = 0.5; %fréquence de coupure
index_fc = ceil((fc*N)/fe);
filtre_ps_haut(1:index_fc) = 0;
filtre_ps_haut(N-index_fc+1:N) = 0;
plot(f,filtre_ps_haut,"linewidth",1.5)
grid on
xlabel('f')
ylabel('Amplitude')
title('Filtre passe-haut')

%Filtrage avec filtre passe-haut

frq_ecg1 = filtre_ps_haut.*y;
ecg1 = ifft(frq_ecg1,"symmetric");
plot(t,ecg);
hold on
plot(t,ecg1+3)
hold on
plot(t,ecg-ecg1+1.5)
grid on
xlabel('t')
ylabel('signal')
title('Représentation du signal avant et aprés le filtrage pass-haut')

%Conception du filtre pass_notch

filtre_passNotch = ones(size(ecg));
fc1 = 50;
index_fc1 = ceil((fc1*N)/fe)+1;
filtre_passNotch(index_fc1) = 0;
filtre_passNotch(N-index_fc1+1) = 0;
plot(f,filtre_passNotch,"linewidth",1.5)
grid on
xlabel('f')
ylabel('Amplitude')
title('Filtre pass-Notch')

%Filtrage avec filtre pass-Notch

frq_ecg2 = filtre_passNotch.*fft(ecg1);
ecg2 = ifft(frq_ecg2, "symmetric");
plot(t,ecg1);
hold on
plot(t,ecg2+3);
hold on
plot(t,ecg-ecg2+1.5);
grid on
xlabel('t')
ylabel('siganl')
title('Représentation du signal avant et aprés le filtrage pass-Notch')

%Partie supprimée

subplot(211)
plot(t,ecg-ecg1)
grid on
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre passe-haut');
subplot(212)
plot(t,ecg-ecg2)
xlim([0 2])
grid on
xlabel('t');
ylabel('signal')
title('Partie supprimée après application du filtre notch');

%Conception du filtre pass-bas

filtre_ps_bas = zeros(size(ecg));
% fc2 = 10;
% fc2 = 20;
fc2 = 30;
index_fc2 = ceil((fc2*N)/fe);
filtre_ps_bas(1:index_fc2) = 1;
filtre_ps_bas(N-index_fc2+1:N) = 1;
plot(f,filtre_ps_bas,"linewidth",1.5)
grid on
xlabel('f'),
ylabel('Amplitude')
title('Filtre pass-bas')

%filtrage avec filtre pass-bas

frq_ecg3 = filtre_ps_bas.*fft(ecg2);
ecg3 = ifft(frq_ecg3, "symmetric");
plot(t,ecg2);
hold on
plot(t,ecg3+3);
hold on
plot(t,ecg-ecg3+1.5);
grid on
xlabel('t')
ylabel('siganl')
title('Représentation du signal avant et aprés le filtrage pass-bas')

%Comparaison des signaux apres application des filtres

subplot(221)
plot(t,ecg)
grid on
xlabel('t')
ylabel('signal')
title('Le signal ecg initiale')

subplot(222)
plot(t,ecg1)
grid on
xlabel('t')
ylabel('signal')
title('Le signal ecg après application du filtre pass-haut')

subplot(223)
plot(t,ecg2)
grid on
xlabel('t')
ylabel('signal')
title('Le signal ecg après application du filtre pass-notch')

subplot(224)
plot(t,ecg3)
grid on
xlabel('t')
ylabel('signal')
title('Le signal après application du filtre pass-bas')

%%
%dernier questions


% Charger le signal ECG
load('ecg3.mat');

% Calculer la fonction d'autocorrelation
[acf, lags] = xcorr(ecg3, ecg3);

% Trouver le premier maximum local après le maximum global
[~, maxIndex] = max(acf);
peakIndex = maxIndex;
for i = maxIndex+1:length(acf)
    if acf(i) < acf(i-1)
        peakIndex = i;
        break;
    end
end

% Calculer la fréquence cardiaque en utilisant le temps entre les pics
heartRate = 1 / (lags(peakIndex) / Fs);

% Afficher la fréquence cardiaque trouvée
disp(['Fréquence cardiaque: ' num2str(heartRate) ' BPM']);






