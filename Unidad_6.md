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
 ### Ejercicio 4
 ```
 function circ(r,x,y)
     xarc(x-r, y+r, r*2, r*2, 0, 360*64)
 endfunction

 function gers(A)
     n = size(A,1)
     c = diag(A)
     r = sum(abs(A), 'c') - abs(c)
     Mx = round(max(c + r) + 1)
     mx = round(min(c - r) - 1)
     My = round(max(r) + 1)
     my = -My
     rect = [mx,my,Mx,My]
     plot2d(real(spec(A)),imag(spec(A)),-1,"031","",rect)
     //replot(rect)
     xgrid()
     for i=1:n
         circ(r(i), c(i), 0)
     end
 endfunction

 function circGersValor(A)
     Gers(A)
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
