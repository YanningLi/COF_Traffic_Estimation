% This function generate quadratic matrix for a number of variabels
% This function is for regularization of a sequence of variables.
% example:
% Regularize q1, q2, q3:
% (q1-q2)^2 + (q2-q3)^2 in matrix form:
% M = quad_matrix(3)
% M = [1  -1  0;
%      -1  2  -1;
%      0  -1  1];

% quad_matrix(1) returns [0];
% quad_matrix(0) returns [] and warning;

function [M] = quad_matrix(n)

    if n==0
        sprintf('Warning: Regularization over an empty variable sequence.\n');
        M = [];
        return
    elseif n==1
        M = [0];
        return
    end

    M = diag([1 2*ones(1,n-2) 1]);
    tr = tril(-1*ones(n,n),-1) - tril(-1*ones(n,n),-2) +...
        triu(-1*ones(n,n),1) - triu(-1*ones(n,n),2);
    M = M + tr;

end