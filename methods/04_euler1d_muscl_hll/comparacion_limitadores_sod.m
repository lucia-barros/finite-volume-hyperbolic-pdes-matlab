%% Comparación de limitadores para MUSCL-Hancock-HLL
% Problema del tubo de choque de Sod
%
% Se comparan:
%   1. Minmod
%   2. Van Leer
%   3. Superbee

clear;
clc;
close all;

%% 1. Parámetros generales

gamma  = 1.4;
xMin   = 0.0;
xMax   = 1.0;
N      = 200;
ng     = 2;
CFL    = 0.45;
tFinal = 0.20;

limitadores = {'minmod', 'vanleer', 'superbee'};

%% 2. Estructura para almacenar los resultados

resultados = struct();

%% 3. Ejecutar MUSCL-Hancock con cada limitador

for k = 1:length(limitadores)

    nombre = limitadores{k};

    fprintf('\nEjecutando MUSCL-Hancock con %s...\n', ...
        upper(nombre));

    [x, rho, u, p, E, iter] = solveSodMUSCL( ...
        N, ng, CFL, tFinal, gamma, ...
        xMin, xMax, nombre);

    resultados.(nombre).x    = x;
    resultados.(nombre).rho  = rho;
    resultados.(nombre).u    = u;
    resultados.(nombre).p    = p;
    resultados.(nombre).E    = E;
    resultados.(nombre).iter = iter;

    fprintf('Iteraciones: %d\n', iter);
    fprintf('Densidad mínima: %.6e\n', min(rho));
    fprintf('Presión mínima: %.6e\n', min(p));

end

%% 4. Representación comparativa

figure;

% ================================================================
% Densidad
% ================================================================

subplot(2,2,1);

plot(x, resultados.minmod.rho, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Minmod');

hold on;

plot(x, resultados.vanleer.rho, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Van Leer');

plot(x, resultados.superbee.rho, ...
    ':', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Superbee');

xlabel('x');
ylabel('\rho');
title('Densidad');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Velocidad
% ================================================================

subplot(2,2,2);

plot(x, resultados.minmod.u, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Minmod');

hold on;

plot(x, resultados.vanleer.u, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Van Leer');

plot(x, resultados.superbee.u, ...
    ':', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Superbee');

xlabel('x');
ylabel('u');
title('Velocidad');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Presión
% ================================================================

subplot(2,2,3);

plot(x, resultados.minmod.p, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Minmod');

hold on;

plot(x, resultados.vanleer.p, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Van Leer');

plot(x, resultados.superbee.p, ...
    ':', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Superbee');

xlabel('x');
ylabel('p');
title('Presión');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Energía total
% ================================================================

subplot(2,2,4);

plot(x, resultados.minmod.E, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Minmod');

hold on;

plot(x, resultados.vanleer.E, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Van Leer');

plot(x, resultados.superbee.E, ...
    ':', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Superbee');

xlabel('x');
ylabel('E');
title('Energía total');
legend('Location','best');
grid on;
hold off;

sgtitle(sprintf(...
    ['Comparación de limitadores MUSCL-Hancock-HLL: ', ...
     'N = %d, t = %.2f'], ...
    N, tFinal));

%% 5. Guardar figura y resultados

exportgraphics(gcf, ...
    'comparacion_limitadores_sod_N200.png', ...
    'Resolution', 300);

save('resultados_limitadores_sod.mat', ...
    'resultados', 'N', 'CFL', 'tFinal', 'gamma');

%% 6. Diferencias máximas respecto a Minmod

difRhoVanLeer = max(abs( ...
    resultados.vanleer.rho - resultados.minmod.rho));

difRhoSuperbee = max(abs( ...
    resultados.superbee.rho - resultados.minmod.rho));

difUVanLeer = max(abs( ...
    resultados.vanleer.u - resultados.minmod.u));

difUSuperbee = max(abs( ...
    resultados.superbee.u - resultados.minmod.u));

difPVanLeer = max(abs( ...
    resultados.vanleer.p - resultados.minmod.p));

difPSuperbee = max(abs( ...
    resultados.superbee.p - resultados.minmod.p));

fprintf('\nComparación respecto al limitador Minmod\n');
fprintf('-----------------------------------------\n');

fprintf('Máxima diferencia de densidad, Van Leer: %.6e\n', ...
    difRhoVanLeer);

fprintf('Máxima diferencia de densidad, Superbee: %.6e\n', ...
    difRhoSuperbee);

fprintf('Máxima diferencia de velocidad, Van Leer: %.6e\n', ...
    difUVanLeer);

fprintf('Máxima diferencia de velocidad, Superbee: %.6e\n', ...
    difUSuperbee);

fprintf('Máxima diferencia de presión, Van Leer: %.6e\n', ...
    difPVanLeer);

fprintf('Máxima diferencia de presión, Superbee: %.6e\n', ...
    difPSuperbee);

%% ================================================================
% FUNCIÓN PRINCIPAL DEL SOLVER
% ================================================================

function [x, rhoFinal, uFinal, pFinal, EFinal, iter] = ...
    solveSodMUSCL(N, ng, CFL, tFinal, gamma, ...
                  xMin, xMax, limiterName)

    %% Malla

    dx = (xMax-xMin)/N;

    x = linspace( ...
        xMin + dx/2, ...
        xMax - dx/2, ...
        N);

    Ntot = N + 2*ng;

    iFis = ng+1 : ng+N;

    %% Condición inicial

    W = zeros(3,Ntot);

    for j = 1:N

        i = j + ng;

        if x(j) < 0.5

            W(:,i) = [1.0; 0.0; 1.0];

        else

            W(:,i) = [0.125; 0.0; 0.1];

        end

    end

    U = primToCons(W, gamma);

    %% Bucle temporal

    t = 0;
    iter = 0;

    while t < tFinal

        iter = iter + 1;

        % Condiciones de contorno
        U = applyTransmissiveBC(U, ng, N);

        % Variables primitivas
        W = consToPrim(U, gamma);

        % Paso temporal CFL
        rho = W(1,iFis);
        vel = W(2,iFis);
        p   = W(3,iFis);

        a = sqrt(gamma*p./rho);

        sMax = max(abs(vel)+a);

        dt = CFL*dx/sMax;

        if t+dt > tFinal
            dt = tFinal-t;
        end

        lambda = dt/dx;

        %% Pendientes limitadas

        slope = zeros(size(W));

        for i = 2:Ntot-1

            deltaLeft = W(:,i)-W(:,i-1);
            deltaRight = W(:,i+1)-W(:,i);

            slope(:,i) = limitedSlope( ...
                deltaLeft, deltaRight, limiterName);

        end

        %% Reconstrucción MUSCL

        WL = W-0.5*slope;
        WR = W+0.5*slope;

        % Protección de positividad
        WL(1,:) = max(WL(1,:),1e-10);
        WR(1,:) = max(WR(1,:),1e-10);

        WL(3,:) = max(WL(3,:),1e-10);
        WR(3,:) = max(WR(3,:),1e-10);

        UL = primToCons(WL,gamma);
        UR = primToCons(WR,gamma);

        %% Predictor Hancock

        FL = physicalFlux(UL,gamma);
        FR = physicalFlux(UR,gamma);

        correction = 0.5*lambda*(FR-FL);

        ULhalf = UL-correction;
        URhalf = UR-correction;

        %% Flujos HLL

        Fnum = zeros(3,Ntot-1);

        for i = 1:Ntot-1

            ULeftInterface  = URhalf(:,i);
            URightInterface = ULhalf(:,i+1);

            Fnum(:,i) = hllFlux( ...
                ULeftInterface, ...
                URightInterface, ...
                gamma);

        end

        %% Actualización conservativa

        Unew = U;

        for i = ng+1:ng+N

            Unew(:,i) = U(:,i) ...
                - lambda*(Fnum(:,i)-Fnum(:,i-1));

        end

        U = Unew;
        t = t+dt;

    end

    %% Solución final

    U = applyTransmissiveBC(U,ng,N);

    Wfinal = consToPrim(U,gamma);

    rhoFinal = Wfinal(1,iFis);
    uFinal   = Wfinal(2,iFis);
    pFinal   = Wfinal(3,iFis);
    EFinal   = U(3,iFis);

end

%% ================================================================
% LIMITADORES
% ================================================================

function slope = limitedSlope(a, b, limiterName)

    slope = zeros(size(a));

    sameSign = (a.*b)>0;

    switch lower(limiterName)

        case 'minmod'

            slope(sameSign) = ...
                sign(a(sameSign)) ...
                .* min(abs(a(sameSign)), ...
                       abs(b(sameSign)));

        case 'vanleer'

            slope(sameSign) = ...
                2*a(sameSign).*b(sameSign) ...
                ./ (a(sameSign)+b(sameSign));

        case 'superbee'

            absA = abs(a(sameSign));
            absB = abs(b(sameSign));

            option1 = min(2*absA,absB);
            option2 = min(absA,2*absB);

            slope(sameSign) = ...
                sign(a(sameSign)) ...
                .* max(option1,option2);

        otherwise

            error('Limitador no reconocido: %s', ...
                limiterName);

    end

end

%% ================================================================
% CONVERSIONES Y FLUJOS
% ================================================================

function U = primToCons(W,gamma)

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = p/(gamma-1)+0.5*rho.*u.^2;

    U = [rho;
         rho.*u;
         E];

end

function W = consToPrim(U,gamma)

    rho = max(U(1,:),1e-12);
    mom = U(2,:);
    E   = U(3,:);

    u = mom./rho;

    p = (gamma-1)*(E-0.5*rho.*u.^2);
    p = max(p,1e-12);

    W = [rho;
         u;
         p];

end

function F = physicalFlux(U,gamma)

    W = consToPrim(U,gamma);

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = U(3,:);

    F = [rho.*u;
         rho.*u.^2+p;
         u.*(E+p)];

end

function Fhll = hllFlux(UL,UR,gamma)

    WL = consToPrim(UL,gamma);
    WR = consToPrim(UR,gamma);

    rhoL = WL(1);
    uL   = WL(2);
    pL   = WL(3);

    rhoR = WR(1);
    uR   = WR(2);
    pR   = WR(3);

    aL = sqrt(gamma*pL/rhoL);
    aR = sqrt(gamma*pR/rhoR);

    SL = min(uL-aL,uR-aR);
    SR = max(uL+aL,uR+aR);

    FL = physicalFlux(UL,gamma);
    FR = physicalFlux(UR,gamma);

    if SL >= 0

        Fhll = FL;

    elseif SR <= 0

        Fhll = FR;

    else

        Fhll = ...
            (SR*FL-SL*FR+SL*SR*(UR-UL)) ...
            /(SR-SL);

    end

end

function U = applyTransmissiveBC(U,ng,N)

    firstPhysical = ng+1;
    lastPhysical  = ng+N;

    U(:,1:ng) = ...
        repmat(U(:,firstPhysical),1,ng);

    U(:,lastPhysical+1:lastPhysical+ng) = ...
        repmat(U(:,lastPhysical),1,ng);

end