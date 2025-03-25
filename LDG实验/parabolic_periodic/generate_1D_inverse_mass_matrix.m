function inv_a = generate_1D_inverse_mass_matrix(mp)
% generate reference mass matrix in[-0.5,0.5] by orthogonal basis polynomial
% inv_a:restore diagonal element
if mp==1
    inv_a=[1,12];
elseif mp==2
    inv_a=[1,12,180];
elseif mp==3
    inv_a=[1,12,180,2800];
end