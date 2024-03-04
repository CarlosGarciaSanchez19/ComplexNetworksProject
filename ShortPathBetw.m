function [d,B] = ShortPathBetw(s,A)
%SHORTPATHBETW is a function that computes the short path betweenness scores
%of all edges in a network starting from a point we will call "s".
% A must be an adjacency matrix representing the network the user wants to
% study. Also this algorithm is thought to be used in an undirected
% network.
% We start by creating two arrays that will represent a weight (important
% part for obtaining the short path betweenness) and a distance that will
% help us during the algorithm.
n = length(A);
w = -ones(1,n);
d = -ones(1,n);
q = -ones(1,n); % queque array
% ALGORITHM of SHORT PATH BETWEENNESS
% 1. Initial vertex is given a distance ds = 0 and a weight ws = 1.
d(s)=0;
w(s)=1;
% 2. Every vertex i adjacent to s is given distance di = ds +1 =1, and weight 
% wi = ws = 1.
i = 1;
j = 2;
l = 1;
k = 0;
q(i) = s;
ts = zeros(1,n);
while i ~= j
    nb = A(q(i),:); % connections of vertex q(i) 
    nbindx = find(nb ~= 0); % neighbor's indices of vertex q(i)
    di = d(q(i));
    w_in = w;
    j_in = j;
    for k=1:length(nbindx)
        if d(nbindx(k)) == -1
            d(nbindx(k)) = di + 1;
            w(nbindx(k)) = A(q(i),nbindx(k))*w(q(i)); % here we multiply by A(i,k) to take into account the case when the adjacency matrix is weighted
            q(j) = nbindx(k);
            j = j+1;
        elseif d(nbindx(k)) == d(q(i)) + 1
            w(nbindx(k)) = w(nbindx(k)) + w(q(i));
        end
    end
    j_out = j;
    w_out = w;
    % Find every "leaf” vertex t, i.e., a vertex such that no paths from s 
    % to other vertices go through t.
    % We do that by saving the position of those vertices from which the path
    % only goes backwards i.e. those for which the weight array stays the same:
    %j == n + 1 &&
    if  sum(w_in == w_out) == sum(w == w)
        ts(l) = q(i);
        l = l + 1;
    end
    if j_in == j_out
        k = i;
    end
    i=i+1;
end
ts = ts(ts ~= 0);
a = ismember(q,ts);
a(1:k-1) = 0;
ts = q(a);
B = zeros(n);
% For each vertex i neighboring t assign a score to the
% edge from t to i of wi/wt.  
i = 1;
while i <= length(ts)
    nbindx = find(A(ts(i),:) ~= 0);
    j = 1;
    while j <= length(nbindx)
        B(ts(i),nbindx(j)) = w(nbindx(j))/w(ts(i));
        j = j + 1;
    end 
    i = i + 1;
end
B = B + transpose(B);
% Now, starting with the edges that are farthest from
% the source vertex s—lower down in a diagram such
% as Fig. 4b—work up towards s. To the edge from
% vertex i to vertex j, with j being farther from s
% than i, assign a score that is 1 plus the sum of
% the scores on the neighboring edges immediately
% below it (i.e., those with which it shares a common
% vertex), all multiplied by wi/wj.
[d_ord,I] = sort(d,'descend');
I2 = I;
for i = 1:length(ts)
    I2(I == ts(i)) = 0;
end
I2 = I2(I2 ~= 0);
i = 1;
while i <= length(I2)
    nbindices = find(A(I2(i),:) ~= 0);
    nbindices_below = nbindices(d(A(I2(i),:) ~= 0) == d_ord(I == I2(i)) + 1);
    nbindices_above = nbindices(d(A(I2(i),:) ~= 0) == d_ord(I == I2(i)) - 1);
    j = 1;
    while j <= length(nbindices_above)
        EB = {(sum(B(I2(i),nbindices_below)) + 1)*w(nbindices_above(j))/w(I2(i)),(sum(B(I2(i),nbindices_below)) + 1)*w(nbindices_above(j))/w(I2(i))};
        [B(I2(i),nbindices_above(j)),B(nbindices_above(j),I2(i))] = EB{:};
        j = j + 1;
    end
    i = i + 1;
end
end

