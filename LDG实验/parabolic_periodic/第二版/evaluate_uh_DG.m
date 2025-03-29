function result = evaluate_uh_DG(uh,x,mp)
% x必须是参考单元上的值

result = 0;
for m = 0:mp

    result = result + uh(m+1)*reference_basis(x,m);

end

end