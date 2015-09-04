% this function takes min of two values which might be empty.
% v = minNonEmpty(v1,v2)

function v = minNonEmpty(v1, v2)

    if isempty(v1) && ~isempty(v2)
        v = min(v2);
    elseif isempty(v2) && ~isempty(v1)
        v = min(v1);
    else
        v = min( min(v1),min(v2) );
    end

end
