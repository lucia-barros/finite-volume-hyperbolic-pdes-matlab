%% Estudio de refinamiento para el tubo de choque de Sod
% Método MUSCL-Hancock-HLL con limitador Van Leer
%
% Mallas consideradas:
% N = 100, 200, 400 y 800
%
% Se representan las soluciones y se calcula una estimación
% de la convergencia mediante diferencias entre mallas consecutivas.

clear;
clc;
close all;

%% 1. Parámetros generales

gamma  = 1.4;

xMin = 0.0;
xMax = 1.0;

ng = 2;

CFL = 0.45;
tFinal = 0.20;

limiterName = 'vanleer';

mallas = [100, 200, 400, 800];

numMallas = length(mallas);

%% 2. Estructura para almacenar los resultados

resultados = struct();

iteraciones = zeros(1,numMallas);
tiemposCPU  = zeros(1,numMallas);

%% 3. Ejecución para cada malla

for k = 1:numMallas

    N = mallas(k);

    fprintf('\n========================================\n');
    fprintf('Ejecutando simulación con N = %d\n',N);
    fprintf('========================================\n');

    tic;

    [x,rho,u,p,E,iter] = solveSodMUSCL( ...
        N,ng,CFL,tFinal,gamma, ...
        xMin,xMax,limiterName);

    tiemposCPU(k) = toc;

    nombreCampo = sprintf('N%d',N);

    resultados.(nombreCampo).x   = x;
    resultados.(nombreCampo).rho = rho;
    resultados.(nombreCampo).u   = u;
    resultados.(nombreCampo).p   = p;
    resultados.(nombreCampo).E   = E;

    iteraciones(k) = iter;

    fprintf('Iteraciones: %d\n',iter);
    fprintf('Tiempo de CPU: %.4f s\n',tiemposCPU(k));
    fprintf('Densidad mínima: %.6e\n',min(rho));
    fprintf('Presión mínima: %.6e\n',min(p));

end

%% 4. Representación de las cuatro mallas

figure;

estilos = {'--','-.',':','-'};
anchuras = [1.2,1.3,1.5,1.8];

% ================================================================
% Densidad
% ================================================================

subplot(2,2,1);
hold on;

for k = 1:numMallas

    N = mallas(k);
    campo = sprintf('N%d',N);

    plot(resultados.(campo).x, ...
         resultados.(campo).rho, ...
         estilos{k}, ...
         'LineWidth',anchuras(k), ...
         'DisplayName',sprintf('N = %d',N));

end

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
hold on;

for k = 1:numMallas

    N = mallas(k);
    campo = sprintf('N%d',N);

    plot(resultados.(campo).x, ...
         resultados.(campo).u, ...
         estilos{k}, ...
         'LineWidth',anchuras(k), ...
         'DisplayName',sprintf('N = %d',N));

end

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
hold on;

for k = 1:numMallas

    N = mallas(k);
    campo = sprintf('N%d',N);

    plot(resultados.(campo).x, ...
         resultados.(campo).p, ...
         estilos{k}, ...
         'LineWidth',anchuras(k), ...
         'DisplayName',sprintf('N = %d',N));

end

xlabel('x');
ylabel('p');
title('Presión');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Energía
% ================================================================

subplot(2,2,4);
hold on;

for k = 1:numMallas

    N = mallas(k);
    campo = sprintf('N%d',N);

    plot(resultados.(campo).x, ...
         resultados.(campo).E, ...
         estilos{k}, ...
         'LineWidth',anchuras(k), ...
         'DisplayName',sprintf('N = %d',N));

end

xlabel('x');
ylabel('E');
title('Energía total');
legend('Location','best');
grid on;
hold off;

sgtitle(sprintf( ...
    ['Estudio de refinamiento: MUSCL-Hancock-HLL ', ...
     'con Van Leer, t = %.2f'], ...
    tFinal));

exportgraphics(gcf, ...
    'refinamiento_sod_muscl_vanleer.png', ...
    'Resolution',300);

%% 5. Cálculo de diferencias entre mallas consecutivas

% Se comparan:
%
% N = 100 con N = 200
% N = 200 con N = 400
% N = 400 con N = 800
%
% La solución fina se restringe a la malla gruesa promediando
% cada pareja de celdas consecutivas.

numComparaciones = numMallas-1;

difRho = zeros(1,numComparaciones);
difU   = zeros(1,numComparaciones);
difP   = zeros(1,numComparaciones);
difE   = zeros(1,numComparaciones);

for k = 1:numComparaciones

    Nc = mallas(k);
    Nf = mallas(k+1);

    campoC = sprintf('N%d',Nc);
    campoF = sprintf('N%d',Nf);

    dxC = (xMax-xMin)/Nc;

    % Solución gruesa
    rhoC = resultados.(campoC).rho;
    uC   = resultados.(campoC).u;
    pC   = resultados.(campoC).p;
    EC   = resultados.(campoC).E;

    % Solución fina
    rhoF = resultados.(campoF).rho;
    uF   = resultados.(campoF).u;
    pF   = resultados.(campoF).p;
    EF   = resultados.(campoF).E;

    % Restricción de la malla fina a la malla gruesa.
    % Para mantener consistencia conservativa, primero se promedian
    % las variables conservadas [rho, rho*u, E] y después se recuperan
    % las variables primitivas. No se promedian u y p directamente.
    UF = [rhoF; rhoF.*uF; EF];

    UFrestringida = 0.5*( ...
        UF(:,1:2:end) + UF(:,2:2:end));

    WFrestringida = consToPrim(UFrestringida, gamma);

    rhoFrestringida = WFrestringida(1,:);
    uFrestringida   = WFrestringida(2,:);
    pFrestringida   = WFrestringida(3,:);
    EFrestringida   = UFrestringida(3,:);

    % Norma L1
    difRho(k) = dxC*sum(abs(rhoC-rhoFrestringida));
    difU(k)   = dxC*sum(abs(uC-uFrestringida));
    difP(k)   = dxC*sum(abs(pC-pFrestringida));
    difE(k)   = dxC*sum(abs(EC-EFrestringida));

end

%% 6. Estimación del orden observado

% Si D_N ~ C*(dx)^p, entonces:
%
% p = log2(D_N/D_2N)

ordenRho = log2(difRho(1:end-1)./difRho(2:end));
ordenU   = log2(difU(1:end-1)./difU(2:end));
ordenP   = log2(difP(1:end-1)./difP(2:end));
ordenE   = log2(difE(1:end-1)./difE(2:end));

%% 7. Mostrar las diferencias

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('DIFERENCIAS L1 ENTRE MALLAS CONSECUTIVAS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf(['Comparación        Densidad        Velocidad       ', ...
         'Presión          Energía\n']);

for k = 1:numComparaciones

    fprintf('%4d - %-4d      %.6e    %.6e    %.6e    %.6e\n', ...
        mallas(k),mallas(k+1), ...
        difRho(k),difU(k),difP(k),difE(k));

end

%% 8. Mostrar los órdenes estimados

fprintf('\n');
fprintf('============================================================\n');
fprintf('ÓRDENES DE CONVERGENCIA ESTIMADOS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf(['Mallas empleadas         Densidad    Velocidad    ', ...
         'Presión     Energía\n']);

for k = 1:length(ordenRho)

    fprintf('%4d-%4d-%4d          %.4f       %.4f       %.4f       %.4f\n', ...
        mallas(k),mallas(k+1),mallas(k+2), ...
        ordenRho(k),ordenU(k),ordenP(k),ordenE(k));

end

%% 9. Tabla de coste computacional

fprintf('\n');
fprintf('============================================================\n');
fprintf('COSTE COMPUTACIONAL\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('N             Iteraciones       Tiempo CPU (s)\n');

for k = 1:numMallas

    fprintf('%-8d      %-10d        %.6f\n', ...
        mallas(k),iteraciones(k),tiemposCPU(k));

end

%% 10. Guardado de resultados

save('resultados_refinamiento_sod.mat', ...
    'resultados','mallas','iteraciones','tiemposCPU', ...
    'difRho','difU','difP','difE', ...
    'ordenRho','ordenU','ordenP','ordenE', ...
    'gamma','CFL','tFinal','limiterName');

%% ================================================================
% SOLVER MUSCL-HANCOCK-HLL
% ================================================================

function [x,rhoFinal,uFinal,pFinal,EFinal,iter] = ...
    solveSodMUSCL(N,ng,CFL,tFinal,gamma, ...
                  xMin,xMax,limiterName)

    %% Malla

    dx = (xMax-xMin)/N;

    x = linspace( ...
        xMin+dx/2, ...
        xMax-dx/2, ...
        N);

    Ntot = N+2*ng;

    iFis = ng+1:ng+N;

    %% Condición inicial

    W = zeros(3,Ntot);

    for j = 1:N

        i = j+ng;

        if x(j)<0.5

            W(:,i) = [1.0;0.0;1.0];

        else

            W(:,i) = [0.125;0.0;0.1];

        end

    end

    U = primToCons(W,gamma);

    %% Bucle temporal

    t = 0;
    iter = 0;

    while t<tFinal

        iter = iter+1;

        % Condiciones de contorno
        U = applyTransmissiveBC(U,ng,N);

        % Variables primitivas
        W = consToPrim(U,gamma);

        % CFL
        rho = W(1,iFis);
        vel = W(2,iFis);
        p   = W(3,iFis);

        a = sqrt(gamma*p./rho);

        sMax = max(abs(vel)+a);

        dt = CFL*dx/sMax;

        if t+dt>tFinal
            dt = tFinal-t;
        end

        lambda = dt/dx;

        %% Pendientes limitadas

        slope = zeros(size(W));

        for i = 2:Ntot-1

            deltaLeft  = W(:,i)-W(:,i-1);
            deltaRight = W(:,i+1)-W(:,i);

            slope(:,i) = limitedSlope( ...
                deltaLeft,deltaRight,limiterName);

        end

        %% Reconstrucción MUSCL

        WL = W-0.5*slope;
        WR = W+0.5*slope;

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
                -lambda*(Fnum(:,i)-Fnum(:,i-1));

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
% LIMITADOR VAN LEER
% ================================================================

function slope = limitedSlope(a,b,limiterName)

    slope = zeros(size(a));

    sameSign = (a.*b)>0;

    switch lower(limiterName)

        case 'vanleer'

            slope(sameSign) = ...
                2*a(sameSign).*b(sameSign) ...
                ./ (a(sameSign)+b(sameSign));

        case 'minmod'

            slope(sameSign) = ...
                sign(a(sameSign)) ...
                .*min(abs(a(sameSign)),abs(b(sameSign)));

        case 'superbee'

            absA = abs(a(sameSign));
            absB = abs(b(sameSign));

            option1 = min(2*absA,absB);
            option2 = min(absA,2*absB);

            slope(sameSign) = ...
                sign(a(sameSign)) ...
                .*max(option1,option2);

        otherwise

            error('Limitador no reconocido: %s',limiterName);

    end

end

%% ================================================================
% FUNCIONES AUXILIARES
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

    if SL>=0

        Fhll = FL;

    elseif SR<=0

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