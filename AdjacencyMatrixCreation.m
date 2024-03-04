function [A] = AdjacencyMatrixCreation(AbstrS)
% Function that creates an adjacency matrix corresponding to a network of
% words.
% The input must be a matrix with rows representing the abstracts and
% columns representing words alone.
    M = length(AbstrS(:,1));
    A = int8(zeros(M));
    for i = 1:(M - 1)
        a = AbstrS(i, AbstrS(i,:) ~= '0');
        for j = (i + 1):M
            b = AbstrS(j, AbstrS(j,:) ~= '0');
            A(i,j) = int8(length(intersect(a,b)));
        end
    end
    A = A + transpose(A);
end

