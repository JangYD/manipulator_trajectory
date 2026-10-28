function [Q, R] = qr_gram_schmidt(matrix)
    [row,col] = size(matrix);
    Q = zeros(row,col);
    R = zeros(col,col);
    V = matrix;

    for i = 1 : col
        R(i,i) = norm(V(:,i));
        Q(:,i) = V(:,i) / R(i,i);
        
        for j = i+1 : col
            R(i,j) = Q(:,i)'*V(:,j);
            V(:,j) = V(:,j) - R(i,j)*Q(:,i);
        end
    end
end
