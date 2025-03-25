function plot_n_DG(uh,P,T)

ng = length(T);
mp = size(T,1)-1;
% uh = reshape(uh,mp+1,ng);
x_list = zeros(1,ng);
y_list = zeros(1,ng);

for n = 1:ng

    phi_local = uh(:,n);
    mesh_nodes = P(:,T(:,n));
    h_local = max(mesh_nodes)-min(mesh_nodes);
    mid_point = mean(mesh_nodes);
    x = linspace(min(mesh_nodes),max(mesh_nodes)-h_local/10,1);
    x_ref = (x -mid_point)/h_local;

    y = zeros(1,1);

    for m = 0:mp
        
        y = y + phi_local(m+1)*reference_basis(x_ref',m);

    end

    x_list(:,n) = x;
    y_list(:,n) = y;

end

x_list = reshape(x_list,1,1*ng);
y_list = reshape(y_list,1,1*ng);


plot(x_list,y_list,'-o');
hold on
plot(x_list,nd(x_list));