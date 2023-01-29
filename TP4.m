clear all
close all
clc

Te=0.0005;
fe = 1/Te;
t = 0:Te:5;

f1=500;
f2=400;
f3=50;

%Représentation temporelle du signal
N=length(t);
x=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t);
plot(t,x);
grid on
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal');

%Représentation fréquentielle du signal
fshift=(-N/2:N/2-1)*fe/N;
y=fft(x);
plot(fshift,fftshift(abs(y)));
grid on
xlabel('f');
ylabel('Amplitude');
title('Représentation fréquentielle en amplitude du signal');

f = (0:N-1)*(fe/N);
K=1;
w=2*pi*f;
H=(K*1i*w/wc)./(1+1i*w/wc);
wc=50;

semilogx(f,abs(H));
grid on
xlabel('f');
ylabel('Amplitude');
title('Représentation du le module de la transmittance complexe');


wc1=500;
wc2=1000;

H1=(K*1i*w/wc1)./(1+1i*w/wc1);
H2=(K*1i*w/wc2)./(1+1i*w/wc2);

G=20*log(abs(H));
G1=20*log(abs(H1));
G2=20*log(abs(H2));

semilogx(f,G,f,G1,f,G2)
grid on


%Filtrage du signal

filtre1 = H.*y;
filtre2 = H1.*y;
filtre3 = H2.*y;

x1=ifft(filtre1,"symmetric");
x2=ifft(filtre2,"symmetric");
x3=ifft(filtre3,"symmetric");

plot(fshift, fftshift(abs(fft(x3))));

subplot(411)
plot(t,x)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal initial');

subplot(412)
plot(t,x1)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal pour wc = 50');

subplot(413)
plot(t,x2)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal wc = 500');

subplot(414)
plot(t,x3)
xlim([0.5 1]);
xlabel('t');
ylabel('x(t)');
title('Représentation temporelle du signal wc = 1000');

[audio,fs]=audioread("test.wav"); %lecture de l'audio

Ts= 1/fs; % Période d'échantillonage
N=length(audio); % le nombre d'échantillons égal à la taille du vecteur music
t=0:Ts:(N-1)*Ts;
plot(t,audio)
title('representation temporelle d un signal audio ')
xlabel('t')
xlim([1 1.2]);
grid on

fshift = (-N/2:N/2-1)*(fs/N); % le pas de discrétisation : fe/N
transfF=fft(audio); % transformée de fourirer
plot(fshift,fftshift(abs(transfF)));
title('representation fréquentielle du signal audio ')
xlabel('f')
ylabel('Amplitude')

fc=4700; %fréquence de coupre, 4700 car l'atténuation est de 0.7
f=(0:N-1)*(fs/N);
k=1;

%Transmittance complexe
H = k./(1+1i*(f/fc).^1000);

%Conception du filtre
H_filter = [H(1:floor(N/2)), flip(H(1:floor(N/2)))];

%Filtrage
f = transfF(1:end-1).*H_filter;

%Signal filtré
filtered_sign=ifft(f,"symmetric");
plot(fshift(1:end-1), fftshift(abs(fft(filtered_sign))));