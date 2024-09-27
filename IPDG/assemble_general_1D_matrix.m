function result = assemble_general_1D_matrix(M,H,R,L,ng,mp)
% assemble marix for the left term
% noting: this assemble is only useful with no nonlinear term,
result = sparse(ng*(mp+1),ng*(mp+1));

% assemble 2 to N-1 elements
for n = 2:ng-1
    
    % 分别表示三种矩阵所在的单元索引坐标
    indexh = n;
    indexr = n+1;
    indexl = n-1;

    result((indexh-1)*(mp+1)+1:indexh*(mp+1),(indexh-1)*(mp+1)+1:indexh*(mp+1)) = M(:,:,indexh)+H(:,:,indexh);
    result((indexh-1)*(mp+1)+1:indexh*(mp+1),(indexr-1)*(mp+1)+1:indexr*(mp+1)) = R(:,:,indexr);
    result((indexh-1)*(mp+1)+1:indexh*(mp+1),(indexl-1)*(mp+1)+1:indexl*(mp+1)) = L(:,:,indexl);

end

% assemble 1,N elemnt
% 1
result(1:mp+1,1:mp+1) =M(:,:,1)+H(:,:,1);
result(1:mp+1,mp+2:2*(mp+1)) = R(:,:,2);

% N
result((ng-1)*(mp+1)+1:ng*(mp+1),(ng-1)*(mp+1)+1:ng*(mp+1)) = M(:,:,ng)+H(:,:,ng);
result((ng-1)*(mp+1)+1:ng*(mp+1),(ng-2)*(mp+1)+1:(ng-1)*(mp+1)) = L(:,:,ng-1);