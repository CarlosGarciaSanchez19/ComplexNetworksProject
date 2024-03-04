function PlotsQAndCommunities(A,Dendrogram,Q,varargin)
m = ceil(sum(A(A ~= 0)));
l = length(Dendrogram);
if isempty(varargin)
    q = Dendrogram{Q == max(Q)};
else
    q = Dendrogram{varargin{1}};
end
[~,e] = ModComp(q,A);
e = e*m;
comm = diag(e);
edges = triu(e) - diag(comm);
[a,b] = find(edges ~= 0);
weights = edges(edges ~= 0);
G = graph(a,b,weights,length(e));
LWidths = 10*G.Edges.Weight./max(G.Edges.Weight);
f1 = figure('Name','Communities');
f2 = figure('Name','Communities with nodes');
f3 = figure('Name','Modularity');
cmap = parula(length(e));
figure(f1)
p = plot(G,'LineWidth',LWidths,'EdgeColor','k');
hold on
title('Plot of communities')
if length(e) <= 10
legend('','Location', 'northwest')
end
Ncolor = 1:length(e);
p.NodeCData = Ncolor;
sz = 50;
Vsize = sz.*comm./max(comm);
Vsize(Vsize == 0) = sz/max(comm);
for i = 1:length(Vsize)
    if order(Vsize(i)) < 0
    Vsize(i) = Vsize(i)*10^(abs(order(Vsize(i))));
    end
end
p.MarkerSize = Vsize;
p.NodeLabel = {};
for i = 1:length(Ncolor)
    plot([NaN NaN],[NaN NaN],'Color',cmap(Ncolor(i),:),'DisplayName',string(Ncolor(i)))
end
figure(f2)
G2 = graph(double(A));
p2 = plot(G2,'Layout','force','EdgeColor','k');
hold on
title('Full Network with nodes colored')
if length(e) <= 10
legend('','Location', 'northwest')
end
Ncolor2 = Ncolor(q);
p2.NodeCData = Ncolor2;
for i = 1:length(Ncolor)
    plot([NaN NaN],[NaN NaN],'Color',cmap(Ncolor(i),:),'DisplayName',string(Ncolor(i)));
end
figure(f3)
plot(1:l,Q,'k-')
title('Modularity of the divisions presented in the dendrogram')
xlabel('Number of divisions (growing number of communities)')
ylabel('Modularity Q')
hold on
plot(1:l,max(Q)*ones(1,l),'b')
legend('Q(l)','Maximum value of Q')
end