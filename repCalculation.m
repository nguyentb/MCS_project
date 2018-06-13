function [ rep_pos, rep_neg ] = repCalculation(rep_pos, rep_neg, exp_mat)
% Reputation Calculation
% Input: current Reputation vector & Experience matrix
% Output: updated Reputation vector

    % Exp_mat matrix will be extracted to attain Exp_mat_pos and Exp_mat_neg
    % Ext_mat(i,j) >= threshold => Ext_mat_pos
    % 0< Ext_mat(i,j) < threshold => Ext_mat_neg
    % Let define threshold

    [~,n] = size(exp_mat);
    threshold = 0.5;
    exp_mat_pos = zeros(n,n);
    exp_mat_neg = zeros(n,n);
    for i=1:n
        for j=1:n
            if (exp_mat(i,j) >= threshold)
                exp_mat_pos(i,j) = exp_mat(i,j);
            elseif (0 < exp_mat(i,j))&&(exp_mat(i,j) < threshold)
                exp_mat_neg(i,j) = 1 - exp_mat(i,j);
            end
        end
    end

%     exp_mat_pos = exp_mat;
%     exp_mat_pos(exp_mat_pos < threshold) = 0;
%     
%     exp_mat_neg = exp_mat;    
%     exp_mat_neg(exp_mat_neg >= threshold) = 0;
%       
    
    % Generate a transition matrix A
    % A(i, j) = 0 if there is no Experience from j to i
    % A(i, j) = exp(j, i)/C(j) where exp(j, i) is experience from j to i;
    % and C(j) is total experience of node j (made by node j to other nodes)

%     c_pos = sum(exp_mat_pos, 2);
%     c_neg = sum(exp_mat_neg, 2);

    c_pos = zeros(1, n);
    c_neg = zeros(1, n);
    for i=1:n
        for j=1:n
            c_pos(1, i) = c_pos(1, i) + exp_mat_pos(i, j);
            c_neg(1, i) = c_neg(1, i) + exp_mat_neg(i, j);
        end
    end
    
    damping = 0.85;
    alpha = (1-damping)/n;
    
    %A as the transition matrix
    A_pos = zeros(n, n);
    A_neg = zeros(n, n);
%     A_pos = (exp_mat_pos.')./fliplr(c_pos);
%     A_pos(isnan(A_pos)) = 0;
%     A_pos(isinf(A_pos)) = 0;
%     
%     A_neg = (exp_mat_neg.')./fliplr(c_neg);
%     A_neg(isnan(A_neg)) = 0;
%     A_neg(isinf(A_neg)) = 0;
    
    for i=1:n
        for j=1:n
            if exp_mat_pos(j, i) > 0
                A_pos(i, j) = exp_mat_pos(j, i)/c_pos(1, j);        
            elseif exp_mat_neg(j, i) > 0
                A_neg(i, j) = exp_mat_neg(j,i)/c_neg(1, j);
            end
        end
    end

    err = 1;
    tol = 1e-10;
    I = ones(n, 1);
    iter = 1;
    while(err>tol)    
        S_pos = damping*A_pos*rep_pos + alpha*I;
        S_neg = damping*A_neg*rep_neg + alpha*I;
        err = norm((S_pos-rep_pos) + norm(S_neg-rep_neg));
        rep_pos = S_pos;
        rep_neg = S_neg;
        iter = iter + 1;
    end
end




