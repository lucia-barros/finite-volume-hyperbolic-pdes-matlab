%% MUSCL-Hancock para las ecuaciones de Euler 1D
% Problema del tubo de choque de Sod
%
% Variables conservadas:
% U = [rho; rho*u; E]
%
% Variables primitivas:
% W = [rho; u; p]

clear;
clc;
close all;

%% 1. Parámetros físicos y numéricos

gamma = 1.4;

xMin = 0.0;
xMax = 1.0;

N = 200;           % Número de celdas físicas
ng = 2;            % Número de celdas fantasma por cada lado

CFL = 0.45;
tFinal = 0.20;

dx = (xMax - xMin)/N;

% Centros de las celdas físicas
x = linspace(xMin + dx/2, xMax - dx/2, N);

% Número total de celdas incluyendo las fantasma
Ntot = N + 2*ng;

%% 2. Inicialización de las variables primitivas

W = zeros(3, Ntot);

% Índices de las celdas físicas
iFis = ng + 1 : ng + N;

for j = 1:N

    i = j + ng;

    if x(j) < 0.5

        % Estado izquierdo del problema de Sod
        rho = 1.0;
        u   = 0.0;
        p   = 1.0;

    else

        % Estado derecho del problema de Sod
        rho = 0.125;
        u   = 0.0;
        p   = 0.1;

    end

    W(:,i) = [rho; u; p];

end

%% 3. Conversión a variables conservadas

U = primToCons(W, gamma);

%% 4. Bucle temporal

t = 0;
iter = 0;

while t < tFinal

    iter = iter + 1;

    %% 4.1 Condiciones de contorno transmisivas

    U = applyTransmissiveBC(U, ng, N);

    %% 4.2 Conversión a primitivas

    W = consToPrim(U, gamma);

    %% 4.3 Cálculo del paso temporal mediante CFL

    rho = W(1,iFis);
    vel = W(2,iFis);
    p   = W(3,iFis);

    a = sqrt(gamma*p./rho);

    sMax = max(abs(vel) + a);

    dt = CFL*dx/sMax;

    % Evitar superar el tiempo final
    if t + dt > tFinal
        dt = tFinal - t;
    end

    lambda = dt/dx;

    %% 4.4 Cálculo de pendientes limitadas

    slope = zeros(size(W));

    for i = 2:Ntot-1

        deltaLeft  = W(:,i)   - W(:,i-1);
        deltaRight = W(:,i+1) - W(:,i);

        slope(:,i) = minmod(deltaLeft, deltaRight);

    end

    %% 4.5 Reconstrucción MUSCL

    WL = W - 0.5*slope;
    WR = W + 0.5*slope;

    % Protección básica frente a valores no físicos
    WL(1,:) = max(WL(1,:), 1e-10);
    WR(1,:) = max(WR(1,:), 1e-10);

    WL(3,:) = max(WL(3,:), 1e-10);
    WR(3,:) = max(WR(3,:), 1e-10);

    %% 4.6 Conversión de estados reconstruidos a conservados

    UL = primToCons(WL, gamma);
    UR = primToCons(WR, gamma);

    %% 4.7 Predictor de Hancock

    FL = physicalFlux(UL, gamma);
    FR = physicalFlux(UR, gamma);

    correction = 0.5*lambda*(FR - FL);

    ULhalf = UL - correction;
    URhalf = UR - correction;

    %% 4.8 Cálculo de flujos HLL en las interfaces

    Fnum = zeros(3, Ntot-1);

    for i = 1:Ntot-1

        % Estado que llega desde la celda izquierda
        ULeftInterface = URhalf(:,i);

        % Estado que llega desde la celda derecha
        URightInterface = ULhalf(:,i+1);

        Fnum(:,i) = hllFlux(ULeftInterface, ...
                            URightInterface, gamma);

    end

    %% 4.9 Actualización conservativa

    Unew = U;

    for i = ng+1:ng+N

        Unew(:,i) = U(:,i) ...
                  - lambda*(Fnum(:,i) - Fnum(:,i-1));

    end

    U = Unew;

    %% 4.10 Actualización del tiempo

    t = t + dt;

    if mod(iter,20) == 0 || abs(t-tFinal) < 1e-12

        fprintf(['Iteración: %5d | Tiempo: %.6f | ', ...
                 'dt: %.3e | Smax: %.4f\n'], ...
                 iter, t, dt, sMax);

    end

end

%% 5. Recuperación de la solución final

U = applyTransmissiveBC(U, ng, N);
Wfinal = consToPrim(U, gamma);

rhoFinal = Wfinal(1,iFis);
uFinal   = Wfinal(2,iFis);
pFinal   = Wfinal(3,iFis);

EFinal = U(3,iFis);

%% 6. Representación gráfica

figure;

subplot(2,2,1);
plot(x, rhoFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('\rho');
title('Densidad');
grid on;

subplot(2,2,2);
plot(x, uFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('Velocidad');
grid on;

subplot(2,2,3);
plot(x, pFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('p');
title('Presión');
grid on;

subplot(2,2,4);
plot(x, EFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('E');
title('Energía total');
grid on;

sgtitle(sprintf(...
    'Tubo de choque de Sod: MUSCL-Hancock-HLL, N = %d, t = %.2f', ...
    N, tFinal));

exportgraphics(gcf, ...
    'sod_muscl_hancock_hll_N200.png', ...
    'Resolution', 300);

%% 7. Comprobaciones básicas

fprintf('\nSimulación terminada correctamente.\n');
fprintf('Número total de iteraciones: %d\n', iter);
fprintf('Densidad mínima: %.6e\n', min(rhoFinal));
fprintf('Presión mínima: %.6e\n', min(pFinal));

%% Guardado de los resultados MUSCL-Hancock

save('resultado_muscl_hancock.mat', ...
    'x', 'rhoFinal', 'uFinal', 'pFinal', ...
    'EFinal', 'N', 'tFinal');

%% ================================================================
% FUNCIONES LOCALES
% ================================================================

function U = primToCons(W, gamma)
% Convierte variables primitivas W=[rho;u;p]
% en variables conservadas U=[rho;rho*u;E].

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = p/(gamma-1) + 0.5*rho.*u.^2;

    U = [rho;
         rho.*u;
         E];

end

function W = consToPrim(U, gamma)
% Convierte variables conservadas U=[rho;rho*u;E]
% en variables primitivas W=[rho;u;p].

    rho = U(1,:);
    mom = U(2,:);
    E   = U(3,:);

    rho = max(rho, 1e-12);

    u = mom./rho;

    p = (gamma-1)*(E - 0.5*rho.*u.^2);
    p = max(p, 1e-12);

    W = [rho;
         u;
         p];

end

function F = physicalFlux(U, gamma)
% Calcula el flujo físico de las ecuaciones de Euler 1D.

    W = consToPrim(U, gamma);

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = U(3,:);

    F = [rho.*u;
         rho.*u.^2 + p;
         u.*(E + p)];

end

function Fhll = hllFlux(UL, UR, gamma)
% Flujo numérico HLL entre dos estados conservados.

    WL = consToPrim(UL, gamma);
    WR = consToPrim(UR, gamma);

    rhoL = WL(1);
    uL   = WL(2);
    pL   = WL(3);

    rhoR = WR(1);
    uR   = WR(2);
    pR   = WR(3);

    aL = sqrt(gamma*pL/rhoL);
    aR = sqrt(gamma*pR/rhoR);

    SL = min(uL-aL, uR-aR);
    SR = max(uL+aL, uR+aR);

    FL = physicalFlux(UL, gamma);
    FR = physicalFlux(UR, gamma);

    if SL >= 0

        Fhll = FL;

    elseif SR <= 0

        Fhll = FR;

    else

        Fhll = ...
            (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL);

    end

end

function value = minmod(a, b)
% Limitador Minmod aplicado componente a componente.

    value = zeros(size(a));

    sameSign = (a.*b) > 0;

    value(sameSign) = sign(a(sameSign)) ...
                    .* min(abs(a(sameSign)), abs(b(sameSign)));

end

function U = applyTransmissiveBC(U, ng, N)
% Aplica dos celdas fantasma transmisivas en cada extremo.

    firstPhysical = ng + 1;
    lastPhysical  = ng + N;

    % Borde izquierdo
    U(:,1:ng) = repmat(U(:,firstPhysical), 1, ng);

    % Borde derecho
    U(:,lastPhysical+1:lastPhysical+ng) = ...
        repmat(U(:,lastPhysical), 1, ng);

end