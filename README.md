# Métodos Numéricos
## Unidad 2
### Ejercicio 1
```
function y=raices(p)
    c = coeff(p,0)
    b = coeff(p,1)
    a = coeff(p,2)
    d = b**2 - 4*a*c
    if(b > 0) then
        y(1) = (-b - sqrt(d))/2*a
        y(2) = (2*c)/(-b-sqrt(d))
    else if(b < 0) then
            y(1) = (2*c)/-b+sqrt(d)
            y(2) = (-b+sqrt(d))/2*a
        end
    end
endfunction
```
### Ejercicio 3-b
```
function z=hornet(x,p)
    if  degree(p) <= 0 then
        z = 0
    else
        ai = coeff(p)
        n = degree(p)
        b(n+1) = ai(n+1)
        for i = n : -1 : 1
            b(i) = ai(i) + x*b(i+1);  
        end
        z = b(1)
    end
endfunction
```
## Unidad 3
### Método Newton
```
function salida = newton(fun, x0, tol, iter)
    deff("y=f(x)", "y=" + fun)
    i = 0
    x1 = x0 - f(x0)/numderivative(f, x0)
    while abs(x1 - x0) > tol && i < iter
        i = i+1
        x0 = x1
        x1 =  x0 - f(x0)/numderivative(f, x0)
    end
    if (abs(x1 - x0) > tol) then 
        disp('Se alcanzó el máximo de iteraciones')
    end
    disp(i)
    salida = x1
endfunction
```
