function R = gmchol(A)
    % Gill-Murray generalized Cholesky decomposition
    % Input matrix A must be non-singular and symmetric
    % Reference: Gill & King (2004)
    
    n = size(A, 1);
    R = eye(n);
    E = zeros(n, n);
    
    norm_A = max(sum(abs(A), 2));
    gamm = max(abs(diag(A)));
    delta = max(eps * norm_A, eps);
    
    for j = 1:n
        theta_j = 0;
        
        for i = 1:n
            sum_val = 0;
            for k = 1:(i-1)
                sum_val = sum_val + R(k, i) * R(k, j);
            end
            R(i, j) = (A(i, j) - sum_val) / R(i, i);
            
            if (A(i, j) - sum_val) > theta_j
                theta_j = A(i, j) - sum_val;
            end
            
            if i > j
                R(i, j) = 0;
            end
        end
        
        sum_val = 0;
        for k = 1:(j-1)
            sum_val = sum_val + R(k, j)^2;
        end
        
        phi_j = A(j, j) - sum_val;
        
        if (j + 1) <= n
            xi_j = max(abs(A(j+1:n, j)));
        else
            xi_j = max(abs(A(n, j)));
        end
        
        beta_j = sqrt(max([gamm, (xi_j / n), eps]));
        
        if delta >= max([abs(phi_j), (theta_j^2) / (beta_j^2)])
            E(j, j) = delta - phi_j;
        elseif abs(phi_j) >= max([delta^2 / (beta_j^2), delta])
            E(j, j) = abs(phi_j) - phi_j;
        elseif (theta_j^2) / (beta_j^2) >= max([delta, abs(phi_j)])
            E(j, j) = (theta_j^2) / (beta_j^2) - phi_j;
        end
        
        R(j, j) = sqrt(A(j, j) - sum_val + E(j, j));
    end
    
    R = R' * R;
end
