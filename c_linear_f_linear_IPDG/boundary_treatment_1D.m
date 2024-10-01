function [Matrix_M,vector_b] = boundary_treatment_1D(boundary_nodes,Matrix_M,vector_b,ng,mp)
% boundary treatment file
% involing periodic,dirichlet,newman,robin boundary condition
% boundary_nodes(1,:) boundary points
% boundary_nodes(2,:) condition class 1==periodic 2==dirichlet 3==newman
% 4==robin
% boundary_nodes(3,:) the elememnts

num_nodes = length(boundary_nodes);

for k = 1:num_nodes

    if boundary_nodes(1,k) == 1
        
        % perodic boundary: u(0,t) = u(2*pi,t),so u_{1/2}^{-}=u_{ng+1/2}*^{-}
        % u_{ng+1/2}^{+}=u_{1}*^{+}

        % 先不考虑不均匀剖分，如果是不均匀剖分，生成矩阵需要现算，会让这个函数的输入参数增多不少
        Matrix_M(1:mp+1,(ng-1)*(mp+1)+1:ng*(mp+1)) = Matrix_M((mp+1)+1:2*(mp+1),1:mp+1);
        Matrix_M((ng-1)*(mp+1)+1:ng*(mp+1),1:mp+1) = Matrix_M(1:mp+1,(mp+1)+1:2*(mp+1));

        vector_b = vector_b;
        break;

    end
end