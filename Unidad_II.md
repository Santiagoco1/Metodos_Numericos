# Unidad 2
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
