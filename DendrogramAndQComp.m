function [Dendrogram,Qmax,MaxQdiv,Q] = DendrogramAndQComp(A)
n = length(A);
m = length(A(A ~= 0));
Adiv = A;
Dendrogram = cell(1,m); % for each index of "Dendrogram" one has a division
% of the network into communities.
[Dendrogram{1,:}] = deal(-ones(n,1));
[Dendrogram{1,1}] = deal(ones(n,1));
i = 1;
l = 2; % index representig the distinct divisions into communities 
% ("Dendrogram" index)
k = 1; % number of communities, at the end k will be equal to n.
while ~isempty(Adiv(Adiv ~= 0))
    %ebsum = zeros(m,1);
    Bsum = zeros(n);
    for j = 1:n
        [~,B] = ShortPathBetw(j,Adiv,true);
        %ebsum = ebsum + eb;
        Bsum = Bsum + B;
    end
    j = 1;
    q = -ones(n,1);
    k_in = k;
    k = 1;
    while sum(q == -1) ~= 0
        d = ShortPathBetw(j,Adiv);
        commindices = d ~= -1;
        q(commindices) = k;
        if sum(q == -1) ~= 0
            I_remain = find(q == -1);
            j = I_remain(1);
            k = k + 1;
        end
    end
    k_out = k;
    if k_in ~= k_out
        Dendrogram{l} = q;
        l = l + 1;
    end
    if order(n) > 1
        Bmax = sort(max(Bsum),'descend');
        Bmax = Bmax(Bmax ~= 0);
        r = 10^(ceil(order(length(Bmax))/2));
        Bdelete = zeros(n,"logical");
        for p = 1:int8(round(length(Bmax)/r))
            BdeleteP = Bsum == Bmax(p);
            Bdelete = Bdelete + BdeleteP;
        end
        Bdelete = logical(Bdelete);
        Adiv(Bdelete) = 0;
    else
        Adiv(Bsum == max(max(Bsum))) = 0;
    end
    i = i + 1;
end
Dendrogram = Dendrogram(1:l-1);
Q = zeros(1,l-1);
i = 1;
while i < l
    Q(i) = ModComp(Dendrogram{i},A);
    i = i + 1;
end
Qmax = max(Q);
MaxQdiv = Dendrogram{Q == max(Q)};
end



