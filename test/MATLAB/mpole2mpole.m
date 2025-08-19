%% Parameters & Data
p = 2;
leng = 2.5;
rho = leng*sqrt(3.0);
th = acos(-leng/rho);
ph = atan2(-leng,-leng);

coeffsO = [ 1,  0.326266+0.777895*1i, -0.610745+0.621683*1i;
            0, -0.50993,              -0.288166-0.687056*1i;                
            0,  0.326266-0.777895*1i, -0.451542;
            0,  0,                    -0.288166+0.687056*1i;
            0,  0,                    -0.610745-0.621683*1i ];

%% Forward rotation
coeffsO_rot = zeros(2*p+1,p+1);
for j = 0:p
    R_th = rot_th(th, j)';
    R_ph = rot_ph(-ph, j);
    coeffsO_rot(1:2*j+1,j+1) = R_th * R_ph * coeffsO(1:2*j+1,j+1);
end

%% Translation along z-axis
coeffsM_rot = zeros(2*p+1,p+1);
for j = 0:p
    for k = -j:j
        % [k, j, A(k,j)]

        % for n=0:min(j+k,j-k)
        %     coeffsM_rot(k_idx,j_idx) = coeffsM_rot(k_idx,j_idx) + ...
        %         coeffsO_rot(k_idx,j_idx-n) * A(0,n) * A(k,j-n) / A(k,j) ...
        %         * rho^n;
        % end

        for n=0:j
            for m=max(k+n-j,-n):min(k+j-n,n)
                if (m ~= 0) 
                    continue;
                end
                
                coeffsM_rot(k+j+1,j+1) = coeffsM_rot(k+j+1,j+1) + ...
                    coeffsO_rot(k-m+j-n+1,j-n+1) * ... 
                    1i^(abs(k)-abs(m)-abs(k-m)) * ...
                    A(m,n) * A(k-m,j-n) * rho^n / A(k,j) ;
            end
        end
    end
end
coeffsM_rot

%% Inverse rotation
coeffsM = zeros(2*p+1,p+1);
for j = 0:p
    R_inv_th = rot_th(th, j);
    R_inv_ph = rot_ph(ph, j);
    coeffsM(1:2*j+1,j+1) = R_inv_ph * R_inv_th * coeffsM_rot(1:2*j+1,j+1);
end
coeffsM

%%
% m=-2;
% n=-2;
% s=0;
% 
% coeff = sum_coeff(th,j,m,n,s)
