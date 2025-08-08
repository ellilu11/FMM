function [D] = rot_th(th, j)
    D = zeros(2*j+1);

    for m = -j:j
        m_idx = m+j+1;
        for n = -j:j
            n_idx = n+j+1;
            s_min = max(m-n,0);
            s_max = min(j+m,j-n);
            for s = s_min:s_max
                D(m_idx,n_idx) = D(m_idx,n_idx) + sum_coeff(th,j,m,n,s);
            end 
            D(m_idx,n_idx) = D(m_idx,n_idx) * ...
                sqrt(factorial(j+n)*factorial(j-n)*factorial(j+m)*factorial(j-m));
        end
    end
end

function [coeff] = sum_coeff(th, j, m, n, s)
    a0 = j+m-s;
    a1 = n-m+s;
    a2 = s;
    a3 = j-n-s;
    coeff = (-1)^(n-m+s) * (cos(th/2)^(a0+a3)) * (sin(th/2)^(a1+a2));
    coeff = coeff ./ (factorial(a0) * factorial(a1) * factorial(a2) * factorial(a3));
end