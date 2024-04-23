function [Q,e] = ModComp(q,A)
% MODCOMP is a function that recieves as an input a network division into
% groups (not fixed number of groups) and the corresponding adjacency matrix
% and returns the corresponding modularity value.
e = zeros(max(q));
m = length(A(A ~= 0));
for i = 1:length(A)
    q_nb = q(A(i,:) ~= 0);
    for j = 1:length(q_nb)
        e(q(i),q_nb(j)) = e(q(i),q_nb(j)) + 1;
    end
end
%e = e + transpose(e);
e = e/m;
Q = trace(e) - sum(e^2,'all');
end

