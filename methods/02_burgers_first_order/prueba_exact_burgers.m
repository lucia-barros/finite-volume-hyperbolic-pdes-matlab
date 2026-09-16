clear
clc
close all

% Dominio
a = 0;
b = 5;
m = 200;

x = linspace(a,b,m+1);

% Posición inicial de la discontinuidad
xb = 2.5;

% Tiempo de representación
t = 1;

%% CASO 1: CHOQUE
wl = 0.5;
wr = 0.1;

w_choque = exact_burgers(x,t,xb,wl,wr);

figure
plot(x,w_choque,'LineWidth',1.5)
grid on
xlabel('x')
ylabel('w(x,t)')
title('Solución exacta de Burgers: choque')

%% CASO 2: RAREFACCIÓN
wl = -1;
wr = 2;

w_rarefaccion = exact_burgers(x,t,xb,wl,wr);

figure
plot(x,w_rarefaccion,'LineWidth',1.5)
grid on
xlabel('x')
ylabel('w(x,t)')
title('Solución exacta de Burgers: rarefacción')