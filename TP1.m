close all
clear all
clc

fe = 10000; %Fréquence d’échantillonnage
te = 1/fe; %Période dééchantillonnage
N = 5000; %Nombre d’échantillons
f = (0:N-1)*(fe/N); %fréquence

t = 0:te:(N-1)*te; %Intervalle
x = 1.2*cos(2*pi*440*t+1.2)+ 3*cos(2*pi*550*t) + 0.6*cos(2*pi*2500*t); %signal périodique x(t)
fshift = (-N/2:N/2-1)*(fe/N);
y = fft(x);

figure(1)
subplot(321)
plot(t,x)
title("Représentation temporelle du signal x(t)")
xlabel("t")
ylabel("x(t)")

subplot(322)
plot(fshift,fftshift(abs(y)))
% plot(f,abs(y))
title("Représentation fréquentielle du signal x(t)")
xlabel("f")
ylabel("Amplitude")

%bruit
xnoise = x+2*randn(size(t));
xnoise2 = x+10*randn(size(t));
fshift = (-N/2:N/2-1)*(fe/N);
ynoise = fft(xnoise);
ynoise2 = fft(xnoise2);

subplot(323)
plot(t,xnoise)
title("Représentation temporelle du signal x(t) bruité")
xlabel("t")
ylabel("bruit")

subplot(324)
plot(fshift,fftshift(abs(ynoise)));
xlabel('f');
ylabel('Amplitude');
title("Représentation fréquentielle en amplitude du signal")

subplot(325)
plot(t,xnoise2)
title("Représentation temporelle du signal x(t) bruité")
xlabel("t")
ylabel("bruit")

subplot(326)
plot(fshift,fftshift(abs(ynoise2)));
xlabel('f');
ylabel('Amplitude');
title("Représentation fréquentielle en amplitude du signal")

%filtre
figure(2)
fc = 2500; %fréquence de coupure
pass_bas = zeros(size(x));
index_fc = ceil((fc*N)/fe);
pass_bas(1:index_fc) = 1;
pass_bas(N-index_fc+1:N) = 1;

subplot(211)
plot(f,pass_bas)
xlabel('f')
ylabel('Amplitude')
title('Filtre pass-bas')

%filtrage
z = pass_bas.*y;%filtre en fct de fréquence
sign_filtr = ifft(z,"symmetric");

subplot(212)
plot(fshift,fftshift(abs(fft(sign_filtr))))
title("Représentation fréquentielle en amplitude du sgnal filté")
xlabel("f")
ylabel("Amplitude")

%signal aprés filtrage

figure(3)

subplot(211)
plot(t,x)
title('Représentation temporelle du signal avant filtrage')
xlabel('t')
ylabel('x(t)')

subplot(212)
plot(t,sign_filtr)
title('Représentation temporelle du signal aprés filtrage')
xlabel('t')
ylabel('x(t)')

s='bluewhale.wav';
[son,Fe]=audioread(s);
sound(son,Fe);

son1 = son(2.45e4: 3.10e4);
sound(son,Fe);

plot(son1);
xlabel('t');
ylabel('son');
title('Représentation temporelle du son roqual bleu')
grid on

N = size(son1);

    %transformation de fourier rapide sur le chant 
    fourier = fft(son1); 
    
    %Densité spectrale du Chant
    Densite_spectrale = abs(fourier).^2/N; 
    
    f = (0:floor(N/2))*(Fe/N)/10;
    plot(f,Densite_spectrale(1:floor(N/2)+1));
    legend("Densité spectrale du chant");
    xlabel("Fréquence (Hz)");
    ylabel("Densité spectrale en puissance");



