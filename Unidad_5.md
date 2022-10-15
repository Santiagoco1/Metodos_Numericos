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
