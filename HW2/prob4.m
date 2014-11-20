clear all; clear all; clc;

n = 1000;
A = triu(rand(n));
B = rand(n);
C = zeros(n);
for i = (1:2:n)
    for j = (1:2:n)
        for k = (i:n-1)
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
            C(i,j+1) = C(i,j+1) + A(i,k)*B(k,j+1);
            
            C(i+1,j) = C(i+1,j) + A(i+1,k+1)*B(k+1,j);
            C(i+1,j+1) = C(i+1,j+1) + A(i+1,k+1)*B(k+1,j+1);
        end
        C(i,j) = C(i,j) + A(i,n)*B(n,j);
        C(i,j+1) = C(i,j+1) + A(i,n)*B(n,j+1);
        
    end
end
D = A*B;
D(n,10)
C(n,10)
