function [D] = rot_ph(ph, j)
    D = zeros(2*j+1);
    for m = -j:j
        m_idx = m+j+1;
        D(m_idx,m_idx) = exp(-1i * m * ph);
    end
end