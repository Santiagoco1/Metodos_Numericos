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
### Ejercicio 4
```
function x = determinante(A)
    
[nA,mA] = size(A) 

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
end;

counter = 0;

a = A

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        a(i,k+1:n) = a(i,k+1:n) - a(k,k+1:n) * a(i,k)/a(k,k);
        a(i,1:k) = 0;              // no hace falta para calcular la solución x
    end;
end;

disp(a)

// Calculo del determinante
x = 1
for k=1:n
    x = x * a(k,k) 
end

endfunction
```
### Ejercicio 6
```
function [x,a] = gausselimPP(A,b)

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    if(abs(a(k,k)) < abs(a(k+1,k))) then
        temp = a(k+1,:);
        a(k+1,:) = a(k,:); 
        a(k,:) = temp;
    end
    
    a(k+1,k+1:n+1) = a(k+1,k+1:n+1) - a(k,k+1:n+1) * a(k+1,k)/a(k,k);      
    a(k+1,k) = 0;
end;

// Sustitución regresiva
x = zeros(1,n);
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    x(i) = (a(i,n+1)-(a(i,1:n) * x'))/a(i,i);
end;
endfunction
```
### Ejercicio 7
```
function [L,U,P]= factorizacion_LU(A)
    U = A
    L = eye(A)
    P = eye(A)
    [m, n] = size(A)
    if m <> n then
        error("La matriz no es cuadrada")
        abort
    end
    for k = 1: m - 1
        kpivot = k; amax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:m);
        U(kpivot,k:m) = U(k,k:m);
        U(k,k:m) = temp;
        temp = P(kpivot,:);
        P(kpivot,:) = P(k,:);
        P(k,:) = temp;
        temp = L(kpivot,1:k-1);
        L(kpivot,1:k-1) = L(k,1:k-1);
        L(k,1:k-1) = temp;
        
        for j = k + 1: m
            L(j,k) = U(j, k)/U(k,k)
            U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
        end
    end;
endfunction
```
### Ejercicio 10
```
function [L,U]= doolittle(A)
    [m, n] = size(A)
    
    for i = 1:n
        for k = i:n
            sum = 0
            for j = 1:i-1
                sum = sum + L(i,j)*U(j,k)
            end
            U(i,k) = A(i,k) - sum
        end
        for k = i:n
            if(i == k)
                L(i,i) = 1
            else
                sum = 0
                for j = 1:i-1
                    sum = sum + L(k,j)*U(j,i)
                end
                L(k,i) = (A(k,i) - sum) / U(i,i)
            end
        end
    end
endfunction

function [y,x,L,U] = ecuacion(A,B)
    [nA,mA] = size(A) 
    [nb,mb] = size(B)
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    [L,U] = doolittle(A)
    
    l = [L B];
    x(n,mb) = 0
    for s = 1:mb
        x(s,s) = l(s,n+s)/l(s,s)
        for i = 2:1:n
            sumk = 0
            for k=1:i
                sumk = sumk + l(i,k)*x(k,s)
            end
            x(i,s) = (l(i,n+s)-sumk)/l(i,i)
        end
    end
    
    u = [U x]
    y(n,mb) = 0
    for s= 1:mb
        y(n,s) = u(n,n+s)/u(n,n)
        for i = n-1:-1:1
            sumk = 0
            for k=i+1:n
                sumk = sumk + u(i,k)*y(k,s)
            end
            y(i,s) = (u(i,n+s)-sumk)/u(i,i)
        end
    end
endfunction
```
### Ejercicio 12-a
```
function [U,ind] = cholesky(A)
eps = 1.0e-8

n = size(A,1)
U = zeros(n,n)

t = A(1,1)
if t <= eps then
    printf('Matriz no definida positiva.\n')
    ind = 0
    return
end

U(1,1) = sqrt(t)
for j = 2:n
    U(1,j) = A(1,j)/U(1,1)
end
    
for k = 2:n
    t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(k,k) = sqrt(t)
    for j = k+1:n
        U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
    end
end
ind = 1

endfunction

function [y,x,U] = sistema(A, B)
    [nA,mA] = size(A) 
    [nb,mb] = size(B)
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    n = nA
    [U] = cholesky(A)
    
    l = [U' B];
    x(n,mb) = 0
    for s = 1:mb
        x(s,s) = l(s,n+s)/l(s,s)
        for i = 2:1:n
            sumk = 0
            for k=1:i
                sumk = sumk + l(i,k)*x(k,s)
            end
            x(i,s) = (l(i,n+s)-sumk)/l(i,i)
        end
    end
    
    u = [U x]
    y(n,mb) = 0
    for s= 1:mb
        y(n,s) = u(n,n+s)/u(n,n)
        for i = n-1:-1:1
            sumk = 0
            for k=i+1:n
                sumk = sumk + u(i,k)*y(k,s)
            end
            y(i,s) = (u(i,n+s)-sumk)/u(i,i)
        end
    end

endfunction
```
