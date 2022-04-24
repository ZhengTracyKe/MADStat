%%%%%%%%%%%%%%%% Construct the coauthorship network %%%%%%%%%%%%%%%%%%%%%%%

load('data/AuthorPaperMatrices.mat')

Nauthor = max(AuPap(:,1));
Npaper = max(AuPap(:,2));
X = sparse(AuPap(:,1), AuPap(:,2), 1, Nauthor, Npaper);
A_weighted = X*X';
A = double(A_weighted >0);   %% edge = coauthored at least 1 paper
A = triu(A,1) + tril(A,-1);   %% remove diagonals
clear AuPap;
clear PapPap;
clear A_weighted;
clear X;

[filename, filepath] = uigetfile('data/author_name.txt');
fullname = fullfile(filepath, filename); %%% get the full path (for fopen)
fid=fopen(fullname,'r','n','UTF-8');
data = textscan(fid,'%s','Delimiter',',');
authorNames = data{1};
fclose(fid);
clear data;
clear name;
clear filename;
clear filepath;
clear fullname;
clear fid;
%%
%%%% A: adjacency matrix of the 47311-by-47311 coauthorship network
%%%% authorNames: names of all the 47311 authors




%%%%%%%%%%%%%% Compute SgnQ pvalues for high-degree nodes %%%%%%%%%%%%%%%%%
N = 1000; 
degs = sum(A,2);
[~, sort_id] = sort(degs, 'descend');
nodes_interest = sort_id(1:N);
SgnQpvals = zeros([N,1]);
N_coauthors = zeros([N,1]);

for i = 1:N
    thisnode = nodes_interest(i);
    neighbors = find(A(thisnode,:)>0);
    N_coauthors(i) = length(neighbors);
    A_person = A(union(neighbors,thisnode), union(neighbors,thisnode));
    [SgnQpvals(i),~,~] = SgnQ(A_person);
end

%%
%%% generate the histogram
subplot(1,2,1)
histogram(N_coauthors,'FaceColor', [0.8500 0.3250 0.0980]);
title('#Coauthors','FontSize',15)
subplot(1,2,2)
histogram(SgnQpvals, 16,'FaceColor',[0 0.4470 0.7410]);
title('SgnQ p-values','FontSize',16)

%%

%%% print the p-values for a few high-degree nodes
name = authorNames(nodes_interest(1:50));
pValue = SgnQpvals(1:50);
Ncoauthor = N_coauthors(1:50); 
table(name, Ncoauthor, pValue)
clear name;
clear Ncoauthor;
clear pValue;


%%

save('DiversityResults_coau.mat', 'authorNames', 'nodes_interest',...
             'N_coauthors', 'SgnQpvals');




%%%%%%%%%%%%%% Draw an example personalized coauthorship network %%%%%%%%%%
person = 'Raymond Carroll';

ID = find(strcmp(authorNames, person));
neighbors = find(A(ID,:)>0);
sets = [neighbors, ID];
m = length(sets);
A_person = A(sets,sets);
Degs = sum(A_person,2);  %%% degrees in the penetwork
names = authorNames(sets);
Nau = sum(A(sets,:),2);  %%% #coauthors (degrees in the whole network - 1)

A_person2 = A(neighbors, neighbors);  %%% remove the hub node
bins = conncomp(graph(A_person2));
bins = [bins';Inf];


pVals = zeros(m,1);
for i = 1:m
    tempnode = sets(i);
    tempsets = union(find(A(tempnode,:)>0), tempnode);
    A_temp = A(tempsets, tempsets); 
    [pVals(i),~,~] = SgnQ(A_temp); %%% p-values based on their own personalized networks
end

%%% sort the nodes and print %%%%%%%%%%%%%%%%%%%
[~, ix] = sort(Nau, 'descend');
A_person = A_person(ix,ix);
names = names(ix);
Nau = Nau(ix);
Degs = Degs(ix);
pVals = pVals(ix);
bins = bins(ix);
table(names, Degs, Nau, pVals, bins)



%%% if needed, run community detection on the giant component %%%%%%%%%
K = 2;
giant = find(bins==mode(bins));
templabels = SCORE(A_person(giant,giant) + eye(length(giant)), K);
labels = zeros(m,1);
labels(giant) = templabels;
table(names, Degs, Nau, pVals, labels)





