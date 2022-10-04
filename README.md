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
### Ejercicio 4
```
function y = dfa(f, v, n, h)
    if (n == 1) 
        y = (f(v+h) - f(v))/h
    else
        y = (dfa(f, v+h, n-1, h) - dfa(f, v, n-1, h))/h;
    end
endfunction

function salida = derivar(fun, v, n, h)
	deff("y=f(x)", "y=" + fun)
	if (n == 0)
	    salida = f(v)
	else
	    salida = dfa(f, v, n, h)
	end
endfunction
```
### Ejercicio 5
```
function y = taylor(fun, n, v, x)
	deff("y=f(x)", "y=" + fun)
	y = f(v)
	for i = 1:1:n
		y = y + (derivar(fun, v, i, 0.1)*((x-v)**i))/factorial(i)
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
### Método Bisección
```
function salida = biseccion(fun, a, b, tol)
    deff("y=f(x)", "y=" + fun)
    if f(a)*f(b) < 0 then
        bool = 1
        while(bool)
            c = (a + b) / 2
            disp(b - c)
            if(abs(b - c) > tol)
                if(f(c) == 0)
                    salida = c
                    bool = 0
                else
                    if (f(a) * f(c)) < 0 then
                        b = c
                    else
                        a = c
                    end
                end
            else
                salida = b
                bool = 0
            end
        end
    end
endfunction
```
### Método Secante
```
function salida=secante(fun, x0, x1, tol)
	deff("y=f(x)", "y=" + fun)
	while abs(x0 - x1) > tol
		x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)))
		x0 = x1
		x1 = x2
	end
	salida = x0
endfunction
```
### Método Newton Multivariable
```
function y=fx(x)
    f1 = (x(1)**2) + (x(1)*(x(2)**3))-9
    f2 = (3*(x(1)**2)*x(2))-4-(x(2)**3)
    y = [f1;f2]
endfunction

function res=fy(x0)
    x = x0(1)
    y = x0(2) 
    f1 = 1+(x**2)-y**2+%e**x*cos(y)
    f2 = 2*x*y+%e**x*sin(y)
    res = [f1; f2]
endfunction

function y=newton(x0, f, it)
    x1 = x0 - ((numderivative(f, x0)**-1)*f(x0))
    i = 0
    while(i < it)
        x0 = x1
        x1 = x0 - ((numderivative(f, x0)**-1)*f(x0))
        i = i + 1
    end
    y = x1
endfunction
```
## Unidad 4
### Ejercicio 2
```
function [x,a, y] = gausselim(A,b)

[nA,mA] = size(A) 
[nb,mb] = size(b)
cont = 0
if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;

for k=1:n-1
    for i=k+1:n
        //a(i, k+1:n+1) = a(i,k+1:n+1) - a(k,k+1:n+1)*a(i,k)/a(k,k); ejercicio d)
        for j=k+1:n+1
            cont = cont +1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        //a(i,1:k) = 0 ejercicio 4
        for j=1:k
            cont = cont +1
            a(i,j) = 0;
        end
    end;
end;

// Sustitución regresiva
cont = cont +1
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    sumk = 0
    //sumk = sumk + a(i,i+1:n)*x(i+1:n); ejercicio d)
    for k=i+1:n
        cont = cont +1
        sumk = sumk + a(i,k)*x(k);
    end;
    cont = cont +1
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
y = cont 
endfunction

// Ejemplos de aplicación
A = [3 -2 -1; 6 -2 2; -9 7 1]
b = [0 6 -1]'

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]'

[x,a] = gausselim(A,b)
disp(x)
disp(a)

A2 = [0 2 3; 2 0 3; 8 16 -1]
b2 = [7 13 -3]'

A2 = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b2 = [-8 -20 -2 4]'

[x2,a2] = gausselim(A2,b2)
```
### Ejercicio 3
```
function [x,a] = gausselim(A,B)

[nA,mA] = size(A) 
[nb,mb] = size(B)
if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;
a = [A B]; // Matriz aumentada

// Eliminación progresiva
n = nA;

for k=1:n-1
    for i=k+1:n
        for j=k+1:n+ mb
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        for j=1:k        
            a(i,j) = 0;  
        end              
    end;
end;

// Sustitución regresiva
x(n,mb) = 0
for s= 1:mb
    x(n,s) = a(n,n+s)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k,s);
        end;
        x(i,s) = (a(i,n+s)-sumk)/a(i,i);
    end;
end

endfunction
```
