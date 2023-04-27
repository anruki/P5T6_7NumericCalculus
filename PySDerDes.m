function [yd, ydd] = PySDerDes(x,y)
% Función que calcula la primera y segunda derivada para un conjunto de
% puntos.
% El algoritmo utilizado es el método de aproximación de la derivada por 
% diferencias finitas.
% Inputs:
%   x, y = vectores fila coordenadas de los puntos de datos
% Outputs:
%   yd = vector fila con la primera derivada en los puntos
%   ydd = vector fila con la segunda derivada en los puntos
    n = numel(x);
    if n ~= numel(y)
        disp("Error: Debe introducir el mismo número de coords y que x")
        return
    end
    d_lagrange = @(x, y, xval) ...
        2*xval-x(2)-x(3) / ((x(1)-x(2))*(x(1)-x(3))) * y(1) + ...
        2*xval-x(1)-x(3) / ((x(2)-x(1))*(x(2)-x(3))) * y(2) + ...
        2*xval-x(1)-x(2) / ((x(3)-x(1))*(x(3)-x(2))) * y(3);

    yd = zeros(1, n);
    yd(1) = d_lagrange(x(1:3), y(1:3), x(1));
    for i = 2: n-1
        yd(i) = d_lagrange(x(i:i+2), y(i:i+2), x(i+1));
    end
    yd(n) = d_lagrange(x(n-2:n), y(n-2:n), x(n));

    ydd = zeros(1, n);
end