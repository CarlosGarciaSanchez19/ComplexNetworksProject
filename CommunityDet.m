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
    error('Second input must be logical true or false.')
end
end

