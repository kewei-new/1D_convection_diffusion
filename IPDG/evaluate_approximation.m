function result = evaluate_approximation(x,uh,mp,P,T)
% evaluate the approximation via polymonial space

uh = reshape(uh,length(uh)/(mp+1),mp+1);
% 二维问题需要把abs(P-x)换成相应的距离公式
min_dis = min(abs(P-x));
%T_index:所求坐标x所在的剖分单元
T_index = [];

index = find(abs(P-x)==min_dis);
for s = 1:size(T,1)
    T_index = [T_index,find(T(s,:)==index)];
end

% 如果是在计算端点处的值，验证内外侧是否有很大的跳跃
result_temp = zeros(1,length(T_index));

for n = 1:length(T_index)

    local_points = P(:,T(:,T_index(n)));
    h_local = max(local_points) - min(local_points);
    mid_point = sum(local_points)/length(local_points);
    x_r = (x-mid_point)/h_local;

    for m = 0:mp

        result_temp(n) = result_temp(n) + uh(T_index(n),m+1)*reference_basis(x_r,m,101,0);

    end
end

if max(result_temp)-min(result_temp) < 1e-8 && n>1
    inf = "keep continous!"
    result = mean(result_temp);
elseif n>1
    inf = "jump largely! discontinous!"
    result = mean(result_temp);
else
    result = result_temp;
end

