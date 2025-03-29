
clear;
k=3;
mo=3;

e = zeros(k,3)
for i = 1:k
    error = solve_1D_DDmodel(mo,i);
    e(i,:) = error;
end

for i = 1:k-1
    s(i,:) = log(e(i,:)./e(i+1,:))/log(2);
end

s