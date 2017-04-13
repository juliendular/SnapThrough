function OOB = OOBmatlab(beta, b, x)

OOB = beta .* x .* (x - b) .* (x - 2*b);

end