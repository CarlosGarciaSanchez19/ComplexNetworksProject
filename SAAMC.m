function [A_w,AbstrS,Data,L] = SAAMC(DataCsv)%,NodeNumber)
% SAAMC (Scientific Abstracts Adjacency Matrix Creation) 
% Datacsv input must be the name of a .csv file. This file has
% to be exported from the web page https://www.scopus.com/results/results.uri?sort=plf-f&src=s&st1=black+holes&sid=bd3be185fe64b74913889ca382dc240a&sot=b&sdt=cl&sl=26&s=TITLE-ABS-KEY%28black+holes%29&origin=resultslist&editSaveSearch=&yearFrom=2015&yearTo=2023&sessionSearchId=bd3be185fe64b74913889ca382dc240a&limit=10
% in the "export" option after your search is done.
% NodeNumber corresponds to the approximate number of nodes that the user
% wants for his network (final number will be less or equal than
% NodeNumber).
    % Importing all data:
    Data = readtable(DataCsv);
    % We extract the abstracts of every article.
    Abstr = string(Data.Abstract);
    % Now we are going to mold our data separating every word and making a
    % string matrix with abstracts as rows and words as columns. Also we filter
    % innecesary abstracts.
    N = 400; %N is an approximation of the maximum number of words in an abstract.
    AbstrS = string(zeros(length(Abstr),N)); 
    k = 1;
    K = false(1,length(Abstr));
    while k <= length(Abstr)
        a = string(zeros(1,N));
        M = length(extract(Abstr(k),lettersPattern));
        if N >= M 
            a(1,1:M) = transpose(extract(Abstr(k),lettersPattern));
            AbstrS(k,:) = a;
            if M == 3 % save the position of those articles without abstract.
                K(k) = 1;
            end
        else 
            K(k) = 1; % save the position of those articles with > 400 words.
        end
        k = k + 1;
    end
    AbstrS(K,:) = []; % eliminate articles without abstract and those with > 400 words.
    Data(K,:) = [];
    % We lower every word in order to make the code work later.
    AbstrS = lower(AbstrS);
    ListOfPronCon = lower(importdata("List of pron, conj, adv.txt"));
    % We want to create an adjagency matrix connecting each article (vertexes)
    % if they have similar words (edges).
    % BUT first we need to throw away the basic words such as: "the", "we",
    % "they", "each", etc. that don't tell us any important connection
    % between articles.
    % We achieve this by localising these words and setting them to '0'.
    for i = 1:length(AbstrS(:,1))
        a = ismember(AbstrS(i,:),ListOfPronCon);
        AbstrS(i,a) = '0';
    end
    % Now the thing is that there may be repetitions of words in each abstract and thus we
    % have to remove them to make the code faster as we only want
    % connections based on different words shared between abstracts.
    % We will do this using the matlab function unique()
    M = length(AbstrS(:,1));
    AbstrSp = string(zeros(M,N));
    for i = 1:M
        a = unique(AbstrS(i,:));
        AbstrSp(i,1:length(a)) = a;
    end
    AbstrS = AbstrSp;
    clear AbstrSp
    AbstrS = sort(AbstrS,2,'descend'); % sort the words leaving 0's at the end.
    L = zeros(M,1); % Vector L to save the number of words in each abstract.
    for i = 1:M
        L(i) = length(AbstrS(i,AbstrS(i,:) ~= '0'));
    end
    N = max(L); % this is the maximum number of different words we have in our abstract data base.
    AbstrS = AbstrS(:,1:N);% we have discarded a large number of '0' elements and so the code will be faster.
    H = histogram(L,'FaceColor','b');
    %M = [];
    %i = 2;
    %while isempty(M)
    %    if sum(L >= N/i) <= NodeNumber
    %        M = i;
    %    end
    %    i = i - 0.001;
    %end
    %AbstrS = AbstrS(L >= N/M,:); % here we discard the articles with less words than 
    % 1/M of the maximum number of words.
    %Data = Data(L >= N/M,:);
    AbstrS = AbstrS(L >= 45 & L <= 105,1:105);
    Data = Data(L >= 45 & L <= 105,:);
    % At this point, with all the abstracts filtered, we can relate each
    % article. We will have a weighted adjagency matrix and the value of
    % each edge will represent how strongly related two articles are.
    A_w = AdjacencyMatrixCreation(AbstrS);
    % Now, as we will have coincidence of words in almost all pair of 
    % articles, we select only those connections with more than a third of 
    % the maximum weight in the adjacency matrix.
    %Amax = max(max(A_w));
    %A_Data = double(A_w >= Amax/3);
    %A_w(A_w <= Amax/3) = 0;
    %Amin = min(A_w(A_w~=0));
    %A_w = double(A_w)/double(Amin);
end