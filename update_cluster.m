function [U_new, center, obj_fcn] = update_cluster(datasample, U, cluster_n, expo)
mf = U.^expo;       % MF matrix after exponential modification
center = mf*datasample./((ones(size(datasample, 2), 1)*sum(mf'))'); % new center
%***********************************************************************
out = zeros(size(center, 1), size(datasample, 1));
if size(center, 2) > 1,
    for k = 1:size(center, 1),
	out(k, :) = sqrt(sum(((datasample-ones(size(datasample, 1), 1)*center(k, :)).^2)'));
    end
else	
    for k = 1:size(center, 1),
	out(k, :) = abs(center(k)-datasample)';
    end
end
dist = out;
%***********************************************************************
obj_fcn = sum(sum((dist.^2).*mf));  % objective function
tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));
