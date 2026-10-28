function [Q,R] = qr_rotation(matrix)
    [row,col]= size(matrix);
    R = matrix;
    Q = eye(row);
    tol = 1e-15;
    for i = 1 : min(row,col)-1
        a = R(i:row, i);
        alpha = -sign(a(1))*hypot(a(1),norm(a(2:end)));
        v = a;
        v(1) = v(1) - alpha;
        tau = (alpha - a(1))/alpha;
        v = v/v(1);

        H = eye(row-i+1) - tau*v*v';
        R(i:end,i:end) = H*R(i:end,i:end);
        Q_tmp = eye(row);
        Q_tmp(i:end,i:end) = H;
        Q = Q*Q_tmp;
        R(abs(R(:,i)) < tol,i) = 0;
    end
end


% function [Q, R] = qr_rotation(matrix)
% 
%     % R = matrix;
%     % Q = eye(row);
%     tol = 1e-15;
%     check_upper = triu(matrix,2);
%     check_loewr = tril(matrix,-2);
%     % if all(check_upper(:) < tol) &&all(check_loewr(:) < tol)
%     %     [Q,R] = triple_diagonal(matrix);
%     % else
%     %     [Q,R] = non_triple_diagonal(matrix);
%     % end 
%     [Q,R] = non_triple_diagonal(matrix);
% 
% end
% 
% function [Q,R] = triple_diagonal(matrix)
%     [row,col]= size(matrix);
%     R = matrix;
%     Q = eye(row);
%     tol = 1e-15;
% 
%     for i = 1 : min(row,col)-1
%         x = R(i,i);
%         y = R(i+1,i);
%         [cs,sn] = given_rotation(x,y);
%         G = [cs sn; -sn cs];
%         % old_bk = R(i+1,i)
% 
%         R(i:i+1,:) = G*R(i:i+1,:)
%         % R(:,i:i+1) = R(:,i:i+1)*G
%         Q_tmp = eye(row);
%         Q_tmp(i:i+1,i:i+1) = G;
%         Q = Q*Q_tmp
%         sign_check = Q*R
%         matrix
%         bk = sign_check(i+1,i);
%         origin_bk = matrix(i+1,i);
%         if sign(bk) ~= sign(origin_bk)
%             R(i,:) = -R(i,:)
%             Q(i,:) = -Q(i,:)
%         end
%         R(abs(R(:,i)) < tol,i) = 0;
%     end
% end
% function [Q, R] = qr_householder(matrix)
%     [row,col]= size(matrix);
%     R = matrix;
%     Q = eye(row);
%     tol = 1e-15;
%     for i = 1 : min(row,col)
%         a = R(i:row, i);
% 
%         if norm(a(2:end)) == 0
%             continue
%         end
% 
%         alpha = sign(a(1));
%         if alpha == 0
%             alpha = 1;
%         end
% 
%         v = a;
%         v(abs(v(:)) < tol) = 0;
%         v_tmp = v;
%         v_tmp(v_tmp == 0) = [];
%         if v(1) == 0
%             if v_tmp(1) < 0
%                 v(1) = v(1) - alpha*norm(a);
%             else
%                 v(1) = v(1) + alpha*norm(a);
%             end
%         else
%             v(1) = v(1) + alpha*norm(a);
%         end
%         norm_v = norm(v);
% 
%         u = v/norm_v;
% 
%         H = eye(row-i+1) - 2*u*u';
%         R(i:end,i:end) = H*R(i:end,i:end);
%         Q_tmp = eye(row);
%         Q_tmp(i:end,i:end) = H;
%         Q = Q*Q_tmp;
%     end
% end