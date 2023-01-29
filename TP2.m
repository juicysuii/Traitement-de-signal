close all
clear all
clc

x = 'phrase.wav';
[data,Fs] = audioread(x);

plot(data)
sound(data,Fs)   %signal normal

sound(data,Fs/2) %pour ralentir
sound(data,Fs*3) %por accélérer

riennesertde = data(36175:81567);
stem(riennesertde)
sound(riennesertde,Fs)

courir = data(81567:104608);
stem(courir)
sound(courir,Fs/2)

ilfaut = data(107373:127405);
stem(ilfaut)
sound(ilfaut,Fs)

partirapoint = data(127405:164516);
stem(partirapoint)
sound(partirapoint,Fs)

sonComplet = [riennesertde, courir, ilfaut, partirapoint];
sound(sonComplet,Fs)

fe = 8192;
te = 5/fe;
N = 5000;
t = (0:N-5)*te;

do1 = 5*cos(2*pi*262*t);
sound(do1,fe)
re = 5*cos(2*pi*294*t);
sound(re,fe)
mi = 5*cos(2*pi*330*t);
sound(mi,fe)
fa = 5*cos(2*pi*349*t);
sound(fa,fe)
sol = 5*cos(2*pi*392*t);
sound(sol,fe)
la = 5*cos(2*pi*440*t);
sound(la,fe)
si = 5*cos(2*pi*494*t);
sound(do1,fe)
do2 = 5*cos(2*pi*523*t);
sound(do2,fe)

music = [do1, re, mi, fa, sol, la, si, do2];
sound(music,fe)

f = (0:N-1)*(fe/N);
frr = fft(music);
signalAnalyzer(fftshift(abs(frr)));
