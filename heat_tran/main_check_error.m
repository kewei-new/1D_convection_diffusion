% check error

clear;clc;
k=3;
mp=2;
for i = 1:k
    error = solver_1D_CD_IPDG(i,mp);
    e(i,:) = error;
end

for i = 1:k-1
    s(i,:) = log(e(i,:)./e(i+1,:))/log(2);
end

s