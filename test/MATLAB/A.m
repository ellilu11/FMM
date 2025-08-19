function [a] = A(m, j)
    a = (-1.0)^j / sqrt(factorial(j-m)*factorial(j+m));
end