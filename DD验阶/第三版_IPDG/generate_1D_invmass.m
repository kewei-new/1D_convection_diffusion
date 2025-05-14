function inv_a = generate_1D_invmass(mp)

if mp==1
    inv_a=[1,12];
elseif mp==2
    inv_a=[1,12,180];
elseif mp==3
    inv_a=[1,12,180,2800];
end