function [ index, sv, sf ] = get_submesh( v , f )

index = [];
sv = [];
sf = [];

numOfVertices = size(v, 2);

inv_index = zeros(1, numOfVertices);
inv_order = 1;
for i = 1:numOfVertices
    if norm(v(1:2, i)) <= 1
        index = [index, i];
        inv_index(i) = inv_order;
        inv_order = inv_order + 1;
        sv = [sv v(:, i)];
    end
end

index_set = unique(index);

numOfFaces = size(f, 2);

for i = 1:numOfFaces
    flag = 1;
    for t = 1:3
        if ~ismember( f(t, i), index_set)
            flag = 0;
        end
    end
    
    if flag == 1
        sf = [sf [inv_index(f(1, i)) ;...
            inv_index(f(2, i))  ; ...
            inv_index(f(3, i))   ] ];
    end
end

end
        
            

