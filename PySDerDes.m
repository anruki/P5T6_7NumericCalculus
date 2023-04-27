function [yd, ydd] = PySDerDes(x, y)
% Función que calcula la primera y segunda derivada para un conjunto de
% puntos. El algoritmo utilizado es el método de aproximación de la 
% derivada por diferencias finitas.
% Inputs:
%   x, y = vectores fila coordenadas de los puntos de datos
% Outputs:
%   yd = vector fila con la primera derivada en los puntos
%   ydd = vector fila con la segunda derivada en los puntos
% Notas: 
%   Por simplicidad se asume que el usuario no introduce x repetidas de 
%   forma que a cada valor de x le corresponde un único valor de y.
%   x es la variable independiente respecto de la que se deriva.

    % Se comprueba que la nube de puntos introducida es válida
    n = numel(x);
    if n ~= numel(y)
        disp("Error: Debe introducir el mismo número de coords y que x")
        return
    end
    if n < 4
        disp("Error: Debe introducir al menos 4 puntos")
        return
    end
    
    % Expresión de la primera derivada para el polinomio de Lagrange de
    % grado 2 para 3 puntos evaluada en xval
    d_lagrange = @(x, y, xval) ...
        2*xval-x(2)-x(3) / ((x(1)-x(2))*(x(1)-x(3))) * y(1) + ...
        2*xval-x(1)-x(3) / ((x(2)-x(1))*(x(2)-x(3))) * y(2) + ...
        2*xval-x(1)-x(2) / ((x(3)-x(1))*(x(3)-x(2))) * y(3);

    % Se calcula la primera derivada en cada punto
    yd = zeros(1, n);
    yd(1) = d_lagrange(x(1:3), y(1:3), x(1));
    for i = 1: n-2
        yd(i+1) = d_lagrange(x(i:i+2), y(i:i+2), x(i+1));
    end
    yd(i+2) = d_lagrange(x(i:i+2), y(i:i+2), x(i+2));

    % Expresión de la segunda derivada para el polinomio de Lagrange de
    % grado 3 para 4 puntos evaluada en xval
    dd_lagrange = @(x, y, xval) ...
        6*xval-2*x(2)-2*x(3)-2*x(4) / ((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))) * y(1) + ...
        6*xval-2*x(1)-2*x(3)-2*x(4) / ((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))) * y(2) + ...
        6*xval-2*x(1)-2*x(2)-2*x(4) / ((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))) * y(3) + ...
        6*xval-2*x(1)-2*x(2)-2*x(3) / ((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))) * y(4);

    % Se calcula la segunda derivada en cada punto
    ydd = zeros(1, n);
    ydd(1) = dd_lagrange(x(1:4), y(1:4), x(1));
    for i = 2: n-2
        ydd(i) = d_lagrange(x(i-1:i+2), y(i-1:i+2), x(i));
    end
    ydd(i+1) = dd_lagrange(x(i-1:i+2), y(i-1:i+2), x(i+1));
    ydd(i+2) = dd_lagrange(x(i-1:i+2), y(n-3:i+2), x(i+2));
end