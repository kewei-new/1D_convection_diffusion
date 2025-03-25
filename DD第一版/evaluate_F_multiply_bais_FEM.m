function result = evaluate_F_multiply_bais_FEM(f_fun,uh,Gauss_coefficient,P,T,Tb,basis_derivative,basis_type)

% mp+1 refers to the number of finite element basis functions i.e. Nbl
ng = length(T);
mp = size(Tb,1)-1;
result = zeros(ng*mp+1,1);
GP = Gauss_coefficient(:,1);
GW = Gauss_coefficient(:,2);
Gpn = length(GW);
uh = reshape(uh,mp+1,length(uh)/(mp+1));
for n = 1:ng

    h_local = P(:,T(2,n))-P(:,T(1,n));
    mid_point = mean(P(:,T(:,n)));
    if nargin(f_fun)==1
        % if right side is a function only depend on x
        for i = 1:mp+1
        % convenient basis function 
            for k = 1:Gpn
                result(Tb(i,n),1) = result(Tb(i,n),1)+GW(k)*feval(f_fun,h_local*GP(k)+mid_point)*basis_reference_FEM(GP(k),basis_derivative,i,basis_type)*h_local^(1-basis_derivative);
            end
        end
    else
        uh_local = uh(:,n);
        for i = 1:mp+1
        % convenient basis function 
            for k = 1:Gpn
            % numerical integration
            % Since f_fun is in a local function, the definition field needs to be reduced
            % Since base_fun is defined in the reference unit, there is no need to restore the
            % Since the integration is performed in a local region, multiplying by h_local is done to perform an affine transformation of the integration region
            % 这里f_fun比较特殊，需要根据具体问题进行修                
                result(Tb(i,n),1) = result(Tb(i,n),1)+GW(k)*feval(f_fun,GP(k),h_local*GP(k)+mid_point,uh_local)*basis_reference_FEM(GP(k),basis_derivative,i,basis_type)*h_local^(1-basis_derivative);
            end
        end        
    end
end
