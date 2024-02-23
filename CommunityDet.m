function [Dendrogram,MaxQdiv,Qmax,Q] = CommunityDet(A,varargin)
%COMMUNITYDET is a function that returns a cell containing the different
%divisons of a network in a dendrogram (Dendrogram), a vector assigning a 
%different number to the vertices that belong to different communities 
%corresponding to the maximum modularity division of the network (MaxQdiv)
%and the value of the maximum modularity (Qmax). Its inputs are A, the
%adjacency matrix, and true or false if the user wants to plot the
%modularity of the different divisions, a plot of the gigant element of the
%network with the vertices colored depending on the community they belong
%to and a plot of the communities and its connections. The division will be
%the one corresponding to the maximum modularity configuration.
if islogical(cell2mat(varargin)) || isempty(varargin)
    switch nargin
        case 1
            [Dendrogram,Qmax,MaxQdiv,Q] = DendrogramAndQComp(A);
        case 2
            [Dendrogram,Qmax,MaxQdiv,Q] = DendrogramAndQComp(A);
            if varargin{1} == true
                PlotsQAndCommunities(A,Dendrogram,Q)
            end
    end
elseif ~isempty(varargin)
    error('Second input musbe logical true or false.')
end

% n = length(A);
% m = length(A(A ~= 0));
% Adiv = A;
% Dendrogram = cell(1,m); % for each index of "Dendrogram" one has a division
% % of the network into communities.
% [Dendrogram{1,:}] = deal(-ones(n,1));
% [Dendrogram{1,1}] = deal(ones(n,1));
% i = 1;
% l = 2; % index representig the distinct divisions into communities 
% % ("Dendrogram" index)
% k = 1; % number of communities, at the end k will be equal to n.
% while i <= m
%     %ebsum = zeros(m,1);
%     Bsum = zeros(n);
%     for j = 1:n
%         [~,B] = ShortPathBetw(j,Adiv);
%         %ebsum = ebsum + eb;
%         Bsum = Bsum + B;
%     end
%     j = 1;
%     q = -ones(n,1);
%     k_in = k;
%     k = 1;
%     while sum(q == -1) ~= 0
%         d = ShortPathBetw(j,Adiv);
%         commindices = d ~= -1;
%         q(commindices) = k;
%         if sum(q == -1) ~= 0
%             I_remain = find(q == -1);
%             j = I_remain(1);
%             k = k + 1;
%         end
%     end
%     k_out = k;
%     if k_in ~= k_out
%         Dendrogram{l} = q;
%         l = l + 1;
%     end
%     Adiv(Bsum == max(max(Bsum))) = 0;
%     i = i + 1;
% end
% Dendrogram = Dendrogram(1:l-1);
% % COMPUTATION OF MODULARITY AND PLOT
% Q = zeros(1,l-1);
% i = 1;
% while i < l
%     Q(i) = ModComp(Dendrogram{i},A);
%     i = i + 1;
% end
% Qmax = max(Q);
% MaxQdiv = Dendrogram{Q(Q == max(Q))};
end

