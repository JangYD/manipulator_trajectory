function [L,U,P] = doolittle_LUP_decomposition(matrix)
    [A_row,A_col] = size(matrix);
    
    P = eye(A_row);
    L = eye(A_row);
    U = matrix;
    for i = 1 : A_col-1
        [~,max_idx] = max(abs(U(i:end,i)));
        U_temp = U(i,:);
        U(i,:) = U(i+max_idx - 1,:);
        U(i+max_idx - 1,:) = U_temp;
        
        pivot = U(i,i);
        factor_v = U(:,i)/pivot;

        P_temp = P(i,:);
        P(i,:) = P(i+max_idx-1,:);
        P(i+max_idx-1,:) = P_temp;
        
        for j = i+1 : A_row 
            if abs(factor_v(j)) < 1e-14
                continue
            end
            L(j,i) = U(j,i)/U(i,i);
            Row_pivot = U(i,i:end);
            Row_target = U(j,i:end);
            U(j,i:end) = Row_target - factor_v(j)*Row_pivot;
        end
                    
        if i-1 == 0 
            continue
        end
        
        L_temp = L(i,1:i-1);
        L(i,1:i-1) = L(i+max_idx - 1,1:i-1);
        L(i+max_idx-1,1:i-1) = L_temp;     
    end
    
    LDU_check = norm(P*matrix -L*U,'fro')/norm(matrix,'fro');
    
    if any(abs(LDU_check(:)) > 1e-12)
        error("is not LUP decomposition")
    end
end