function [L,U,x]=Thomas_Algorithm(matrix,b)
    A = matrix;
    %%Tridiagonal_check 
    upper_A = triu(A,2);
    lower_A = tril(A,-2);
    if abs(norm(upper_A)) > 1e-14
        error('A is not Tridiagonal(upper)')
    end
    if abs(norm(lower_A)) > 1e-14
        error('A is not Tridiagonal(lower)')
    end
    [row,~] = size(A);
    L = eye(row);
    U = A;
    U = triu(U,1);
    U(1,1) = A(1,1);
    y = zeros(row,1);
    x = zeros(row,1);
    for i = 1 : row -1
        L(i+1, i) = A(i+1,i)/U(i,i);
        U(i+1,i+1) = A(i+1,i+1) - L(i+1,i)*A(i,i+1);
    end

    for i = 1 : row
        if i == 1 
            y(i) = b(i);
            continue
        end
        y(i) = b(i) - L(i,i-1)*y(i-1);
    end
    check_Ly = L*y;
    if norm(check_Ly - b', inf)  > 1e-14
        error("Ly ~= b")
    end

    for i = row : -1 :  1
        if i == row
            x(i) = y(i)/U(row,row);
            continue
        end
        x(i) = (y(i)-U(i,i+1)*x(i+1))/U(i,i);
    end
    check_LUx = L*U*x;
    if norm(check_LUx - b',inf) > 1e-14
        error("LUx ~= b")
    end
end