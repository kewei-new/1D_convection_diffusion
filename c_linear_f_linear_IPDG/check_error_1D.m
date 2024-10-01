% check_error_1D.m
% check the order of error
clear;clc;
k=3;

for i = 1:k
    error = solver_1D_CD_IPDG(i);
    e(i,:) = error;
end

for i = 1:k-1
    s(i,:) = log(e(i,:)./e(i+1,:))/log(2);
end

s