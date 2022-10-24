function R = buildCovarianceMatrix(rho,samplenum)
c = zeros(1,5);

for i = 1:5
    
    if samplenum < 2000
    sample = rho(i,samplenum-599:samplenum);
    else
     sample = rho(i,samplenum-1499:samplenum);
    end

    c(i) = cov(sample);
end

R = diag(c);
end

