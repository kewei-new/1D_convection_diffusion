function plot_phi_FEM(phi,P,T,Tb,basis_type)

ng = length(T);
mp = size(Tb,1)-1;
phi = phi(Tb);
x_list = zeros(10,ng);
y_list = zeros(10,ng);

for n = 1:ng

    phi_local = phi(:,n);
    mesh_nodes = P(:,T(:,n));
    h_local = max(mesh_nodes)-min(mesh_nodes);
    mid_point = mean(mesh_nodes);
    x = linspace(min(mesh_nodes),max(mesh_nodes)-h_local/10,10);
    x_ref = (x -mid_point)/h_local;

    y = zeros(10,1);

    for m = 1:mp+1
        
        y = y + phi_local(m)*basis_reference_FEM(x_ref',0,m,basis_type);

    end

    x_list(:,n) = x;
    y_list(:,n) = y;

end

x_list = reshape(x_list,1,10*ng);
y_list = reshape(y_list,1,10*ng);

% plot(x_list,y_list,'-o',x_list,temp_exact_fun(x_list),'-*');
plot(x_list,y_list,'-o');