load('data/YearlyCiteeNetworks.mat');


%%% aggregate the citee networks from 1991 to 2000
startYear = 1991;
endYear = 2000;
k_min = find(Year_all==startYear);
k_max = find(Year_all==endYear);
rawAdj = CiteeAdjYearly{k_min};
for k = (k_min+1):k_max
    rawAdj = rawAdj + CiteeAdjYearly{k};
end
clear 'CiteeAdjYearly';
clear 'Year_all';

degrees = sum(rawAdj);
keepNode = (degrees > 60);   %%% node screening
A = double(rawAdj(keepNode,keepNode)>=2);   %%% edge screening
A = triu(A,1)+tril(A, -1);   %%% remove the diagonal
keepNodeID = find(keepNode);

G = graph(A);
bins = conncomp(G);
unique(bins)    %%% it is seen that there is only one giant component

save('CiteeAdjFinal.mat','A','keepNodeID');