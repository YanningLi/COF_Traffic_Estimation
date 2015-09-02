% unique_tol
% this function has the same functionality as unique(), but with a
% tolerance
% also works for matrix

function [uniq] = unique_tol(array,tol)

uniq = unique(array);

len = length(uniq);
index = (uniq(2:len)-uniq(1:len-1))<tol;

% remove the value
if sum(index)~=0
    uniq(index) = [];
end

end