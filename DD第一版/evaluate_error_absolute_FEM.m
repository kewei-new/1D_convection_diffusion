function result = evaluate_error_absolute_FEM(exact_fun,uh,P,T,Tb,ng,mp,basis_type)
% 只计算左端点
uh = uh(Tb);

result = 0;
for n = 1:ng
    
    mesh_nodes = P(:,T(:,n));
    h_local = max(mesh_nodes) - min(mesh_nodes);
    mid_point = sum(mesh_nodes)/length(mesh_nodes);
    x = linspace(min(mesh_nodes),max(mesh_nodes)-h_local/10,6);
    x_ref = (x -mid_point)/h_local;
    uh_local = uh(:,n);
    
    uh_a = 0;
    for m = 1:mp+1
        uh_a = uh_a + uh_local(m)*basis_reference_FEM(x_ref',0,m,basis_type);
    end
    temp = abs(feval(exact_fun,x')-uh_a);

    result = max(max(result,temp));

end