function plot_phi_FEM(phi,P,T,Tb,basis_type,derivative)

ng = length(T);
Tb = Tb(:,2:ng+1);
mp = size(Tb,1)-1;
phi = phi(Tb);
x_list = zeros(1,ng);
y_list = zeros(1,ng);

for n = 1:ng

    phi_local = phi(:,n);
    mesh_nodes = P(:,T(:,n));
    h_local = max(mesh_nodes)-min(mesh_nodes);
    mid_point = mean(mesh_nodes);
    x = linspace(min(mesh_nodes),max(mesh_nodes)-h_local/10,1);
    x_ref = (x -mid_point)/h_local;

    y = zeros(1,1);

    for m = 1:mp+1
        
        y = y + phi_local(m)*reference_basis_FEM(x_ref',derivative,m,basis_type)*h_local^(-derivative);

    end

    x_list(:,n) = x;
    y_list(:,n) = y;

end

x_list = reshape(x_list,1,1*ng);
y_list = reshape(y_list,1,1*ng);

% plot(x_list,y_list,'-o',x_list,temp_exact_fun(x_list),'-*');
plot(x_list,y_list,'-o');