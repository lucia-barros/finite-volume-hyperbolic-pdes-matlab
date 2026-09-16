clear
clc
close all

%% 1. Dominio espacial

a = 0;
b = 5;

m = 100;

x = linspace(a,b,m+1);
dx = x(2)-x(1);

%% 2. Datos del problema de Riemann

xb = 0.5*(a+b);

wl = 0.5;
wr = 0.1;

w0 = zeros(size(x));
w0(x <= xb) = wl;
w0(x > xb)  = wr;

%% 3. Parámetros temporales

CFL = 0.8;

velocidad_maxima = max(abs([wl,wr]));

dt = CFL*dx/velocidad_maxima;

t_final = 2;

% Ajustamos dt para llegar exactamente a t_final
nt = ceil(t_final/dt);
dt = t_final/nt;

lambda = dt/dx;

fprintf('dx = %.6f\n',dx)
fprintf('dt = %.6f\n',dt)
fprintf('Número de pasos temporales = %d\n',nt)
fprintf('CFL real = %.6f\n',velocidad_maxima*dt/dx)

%% 4. Inicialización de todos los métodos

w_godunov = w0;
w_lxf = w0;
w_q = w0;
w_nc = w0;

t = 0;

%% 5. Evolución temporal

for n = 1:nt

    w_godunov = godunov_burgers(w_godunov,lambda);
    w_lxf     = lax_friedrichs_burgers(w_lxf,lambda);
    w_q       = qscheme_burgers(w_q,lambda);
    w_nc      = no_conservativo_burgers(w_nc,lambda);

    t = t+dt;
end

%% 6. Solución exacta

w_exacta = exact_burgers(x,t,xb,wl,wr);

%% 7. Representación de las soluciones

figure

plot(x,w_exacta,'k','LineWidth',2)
hold on

plot(x,w_godunov,'--','LineWidth',1.3)
plot(x,w_lxf,'-.','LineWidth',1.3)
plot(x,w_q,':','LineWidth',1.5)
plot(x,w_nc,'LineWidth',1.3)

grid on
xlabel('x')
ylabel('w(x,t)')
title(sprintf('Ecuación de Burgers, t = %.2f',t))

legend('Solución exacta', ...
       'Godunov', ...
       'Lax-Friedrichs', ...
       'Q-scheme', ...
       'No conservativo', ...
       'Location','best')

%% 8. Errores punto a punto

figure

plot(x,abs(w_exacta-w_godunov),'LineWidth',1.3)
hold on
plot(x,abs(w_exacta-w_lxf),'LineWidth',1.3)
plot(x,abs(w_exacta-w_q),'LineWidth',1.3)
plot(x,abs(w_exacta-w_nc),'LineWidth',1.3)

grid on
xlabel('x')
ylabel('|w_{exacta}-w_{numérica}|')
title(sprintf('Errores puntuales, t = %.2f',t))

legend('Godunov', ...
       'Lax-Friedrichs', ...
       'Q-scheme', ...
       'No conservativo', ...
       'Location','best')

%% 9. Cálculo de errores

error_inf_godunov = max(abs(w_exacta-w_godunov));
error_inf_lxf     = max(abs(w_exacta-w_lxf));
error_inf_q       = max(abs(w_exacta-w_q));
error_inf_nc      = max(abs(w_exacta-w_nc));

error_L1_godunov = dx*sum(abs(w_exacta-w_godunov));
error_L1_lxf     = dx*sum(abs(w_exacta-w_lxf));
error_L1_q       = dx*sum(abs(w_exacta-w_q));
error_L1_nc      = dx*sum(abs(w_exacta-w_nc));

fprintf('\nErrores en t = %.4f\n',t)
fprintf('-------------------------------------------\n')
fprintf('Godunov:\n')
fprintf('  Error infinito = %.6e\n',error_inf_godunov)
fprintf('  Error L1       = %.6e\n',error_L1_godunov)

fprintf('Lax-Friedrichs:\n')
fprintf('  Error infinito = %.6e\n',error_inf_lxf)
fprintf('  Error L1       = %.6e\n',error_L1_lxf)

fprintf('Q-scheme:\n')
fprintf('  Error infinito = %.6e\n',error_inf_q)
fprintf('  Error L1       = %.6e\n',error_L1_q)

fprintf('No conservativo:\n')
fprintf('  Error infinito = %.6e\n',error_inf_nc)
fprintf('  Error L1       = %.6e\n',error_L1_nc)