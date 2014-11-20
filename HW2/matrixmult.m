n = 1024;
A = triu(rand(n));
B = rand(n);
C = zeros(n);
for k = 1:n
    for i = 1:k
        for j = 1:n
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
        end
    end
end
A*B==C