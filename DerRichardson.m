function [dfdx,d2fdx] = DerRichardson(f,x0,hin,N)
% Función que calcula la derivada primera y segunda de una función anónima
% con el método de Richardson
% Entradas:
% f: función dada en forma anónima (f = @(x))
% x0: punto donde se calcula la derivada
% hin: tamaño de paso inicial
% N: orden de extrapolación
%
% Salidas:
% dfdx: matriz de extrapolación para la 1ª derivada
% d2fdx: matriz de extrapolación para la 2ª derivada

% Inicializar las matrices de extrapolación con ceros
dfdx = zeros(N, N+1);
d2fdx = zeros(N, N+1);

% Calcular las primeras columnas de las matrices de extrapolación
for i = 1:N
    h = hin/2^(i-1);
    dfdx(i,1) = h;
    d2fdx(i,1) = h;
    dfdx(i,2) = (feval(f,x0+h) - feval(f,x0-h))/(2*h); % diferencias centrales de dos puntos
    d2fdx(i,2) = (feval(f,x0+h) - 2*feval(f,x0) + feval(f,x0-h))/(h^2); % diferencias centrales de tres puntos
end

% Calcular las extrapolaciones
for j = 3:N+1
    for i = j-1:N
        dfdx(i,j) = dfdx(i,j-1) + (dfdx(i,j-1) - dfdx(i-1,j-1))/(4^(j-2)-1);
        d2fdx(i,j) = d2fdx(i,j-1) + (d2fdx(i,j-1) - d2fdx(i-1,j-1))/(4^(j-2)-1);
    end
end