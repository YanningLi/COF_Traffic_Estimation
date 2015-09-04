% This file is used for nonlinear mapping

% k is the original matrix
% original interval in the form 
% [interval_1_start_point, interval_1_end_point;
%   ...interval_n_start_point, interval_n_end_point];
% target_interval is in the same form. They must have same dimensions.

% Example: 
% k_trans = mapping(k,[0 k_c; k_c k_m],[0 0.5*k_m; 0.5*k_m k_m]);
% mapping (0~k_c)(k_c~k_,m) ==> (0~0.5*k_m)(0.5*k_m~k_m)
% You can define inconinuous target intervals to not use certain portion of
% the colorbar

function [k_trans] = mapping(k,original_interval, target_interval)

k_trans = 0.0*k;
for i = 1:size(original_interval,1)
    isInDomain = (k>=original_interval(i,1) & k<=original_interval(i,2)+10e-6);
    k_trans(isInDomain) = target_interval(i,1) + ...
        (k(isInDomain)-original_interval(i,1))*(target_interval(i,2)-target_interval(i,1))...
        /(original_interval(i,2)-original_interval(i,1));
end
