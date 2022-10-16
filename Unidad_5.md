# Unidad 5
### Método de Jacobi
```
function x = jacobi(A, b, x0, eps)
    n = size(A,1)
    x = x0
    xk = x
    for i = 1:n
        sum = 0
        for j = 1:n
            if (i<>j)
                sum = sum + A(i,j)*xk(j)
            end
        end
        x(i) = 1/A(i,i)*(b(i)-sum)
    end
    while (abs(norm(x-xk)) > eps)
        xk = x
        for i = 1:n
            sum = 0
            for j = 1:n
                if (i<>j)
                    sum = sum + A(i,j)*xk(j)
                end
            end
            x(i) = 1/A(i,i)*(b(i)-sum)
        end
    end
endfunction
```
### Método de Gauss-Seidel
```
function x = gauss_seidel(A, b, x0, eps)
    n = size(A,1)
    x = x0
    xk = x
    for i = 1:n
        sum = 0
        for j = 1:n
            if (i<>j)
                sum = sum + A(i,j)*x(j)
            end
        end
        x(i) = 1/A(i,i)*(b(i)-sum)
    end
    while (abs(norm(x-xk)) > eps)
        xk = x
        for i = 1:n
            sum = 0
            for j = 1:n
                if (i<>j)
                    sum = sum + A(i,j)*x(j)
                end
            end
            x(i) = 1/A(i,i)*(b(i)-sum)
        end
    end
endfunction
```
### Método SOR
```
function w = omega_SOR(A)
    n = size(A,1)
    T_j = eye(n,n) - diag(1./diag(A))*A
    autovalores = spec(T_j)
    rho = max(abs(autovalores))
    w = 2/(1+sqrt(1-rho**2))
endfunction


function [x, iter] = SOR(A, b, x0, tol)
    n = size(A,1)
    w = omega_SOR(A)
    x = x0 + ones(n,1)
    iter = 0
    while (abs(norm(x-x0)) > tol)
        x = x0
        x0(1) = (1-w)*x0(1) + (w/A(1,1)) * (b(1) - A(1,2:n) * x0(2:n))
        for i = 2:n-1
            x0(i) = (1-w) * x0(i) + (w/A(i,i)) * (b(i) - A(i,1:i-1) * x0(1:i-1) - A(i,i+1:n) * x0(i+1:n))
        end
        x0(n) = (1-w) * x0(n) + (w/A(n,n)) * (b(n) - A(n,1:n-1) * x0(1:n-1))
        iter = iter + 1
    end
    x = x0
endfunction
```
