function result = evaluate_uh_DG(x,uh,mp)
% 计算局部区间在参考单元上的u_h
% x:属于[-0.5,0.5]
% uh:基函数系数
% mp:基函数次数(也可以省略，直接由uh生成)
result = 0;

for m = 0:mp

    result = result + uh(m+1)*reference_basis(x,m);

end
