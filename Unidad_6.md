 # Unidad 6
 ### Polinomio Caracteristico
 ```
 function p = pol_car(A)
    x = poly([0 1], "x", "coeff")
    n = size(A, 1)
    I = eye(n,n)
    p = det(A - x * I)
 endfunction
 ```
 ### Método Potencia con Tolerancia
 ```
 function [z,l] = pot(A,z,esp)
    n = size(A,1)
    a = real(spec(A))
    w = A * z
    while(norm(w-z) < esp)
        z = w/norm(w)
        w = A * z
    end
    k = 1
    while(w(k) <> 0) & (k < n)
        k = k+1
    end
    if(w(k) == 0) then
        disp("La matriz es singular")
    else
        l = w(k) / z(k)
    end
    z = w/norm(w)
 endfunction
 ```
 ### Método Potencia con Iteraciones
 ```
 function [d,l] = potb(A,z,iter)
    n = size(A,1)
    a = real(spec(A))
    w = A * z
    for i = 1:iter
        z = w/norm(w)
        w = A * z
    end
    k = 1
    while(w(k) <> 0) & (k < n)
        k = k+1
    end
    if(w(k) == 0) then
        disp("La matriz es singular")
    else
        l = w(k) / z(k)
    end
    z = w/norm(w)
endfunction
 ```
