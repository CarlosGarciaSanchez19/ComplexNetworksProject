load('FinalData2.mat')
clear L
f1 = figure('Name','Histogram1');
f2 = figure('Name','Histogram2');
f3 = figure('Name','Histogram3');
f4 = figure('Name','WeightedNet1');
f5 = figure('Name','WeightedNet2');
figure(f1)
histogram(A_w(A_w ~= 0),'BinEdges',1:104,'FaceColor','r');
title('Histogram of edges from 1 to 104')
xlabel('Edge value')
ylabel('Number of edges')
figure(f2)
histogram(A_w(A_w ~= 0),'BinEdges',50:104,'FaceColor','r');
title('Histogram of edges from 50 to 104')
xlabel('Edge value')
ylabel('Number of edges')
Ap_w = A_w;
Ap_w(Ap_w <= 50) = 0;
A_w(sum(Ap_w(:,:)) == 0,:) = []; A_w(:,sum(Ap_w(:,:)) == 0) = [];
AbstrS(sum(Ap_w(:,:)) == 0,:) = []; Data(sum(Ap_w(:,:)) == 0,:) = [];
figure(f3)
histogram(A_w(A_w ~= 0),'BinEdges',15:104,'FaceColor','g');
title('Histogram from 15 to 104 after removing nodes with connections <= 50')
xlabel('Edge value')
ylabel('Number of edges')
clear Ap_w
A_w(A_w < 18) = 0;
figure(f4)
G2 = graph(double(A_w));
plot(G2,'Layout','force','EdgeColor','k','NodeColor','b');
title('Weighted network after edges with value < 18 removal')
d = ShortPathBetw(148,A_w);
A_w(d == -1,:) = []; A_w(:,d == -1) = [];
AbstrS(d == -1,:) = []; Data(d == -1,:) = [];
figure(f5)
G2 = graph(double(A_w));
plot(G2,'Layout','force','EdgeColor','k','NodeColor','r');
title('Final weighted network')
clear G2 d f1 f2 f3 f4 f5