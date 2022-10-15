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
