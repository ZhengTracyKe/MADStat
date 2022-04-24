%%%%%%%%%%%%%%%% Construct the coauthorship network %%%%%%%%%%%%%%%%%%%%%%%

load('data/AuthorPaperMatrices.mat')

Nauthor = max(AuPap(:,1));
Npaper = max(AuPap(:,2));
X = sparse(AuPap(:,1), AuPap(:,2), 1, Nauthor, Npaper);
A_weighted = X*X';
A = double(A_weighted >=3);   %% edge = coauthored at least 3 papers
A = triu(A,1) + tril(A,-1);   %% remove diagonals
clear AuPap;
clear A_weighted;

G = graph(A);
bins = conncomp(G);
k = mode(bins);  %% find the largest component
keepNodeID = find(bins==k);
A = A(keepNodeID, keepNodeID);
n = size(A,1);


[filename, filepath] = uigetfile('data/author_name.txt');
fullname = fullfile(filepath, filename); %%% get the full path (for fopen)
fid=fopen(fullname,'r','n','UTF-8');
data = textscan(fid,'%s','Delimiter',',');
name = data{1};
fclose(fid);
authorNames = name(keepNodeID);
clear data;
clear name;

save('CoauAdjFinal.mat', 'A', 'authorNames');





%%%%%%%%%%%%%% Apply SCORE for first-layer community detection%%%%%%%%%%%%%
K = 6;
templabels = SCORE(A + eye(n), K);  %%%% the diagonals are added back for numerical 
                                    %%%%stabilization of eigenvector calculation 
summary = tabulate(templabels);

[~, ix] = sort(summary(:,2), 'descend');   %%% permute the order of communities
pemulabel = [5,1,4,3,6,2];
labels_firstlayer = zeros(n,1);
for k =1:6 
    labels_firstlayer(templabels==ix(k)) = pemulabel(k);
end
tabulate(labels_firstlayer)
clusterNames_firstlayer = {'Nonparametric Statistics', 'Biostatistics (Europe)',...
                          'Mathematical Statistics', 'Biostatistics (UNC)', ...
                          'Semi-parametric Statistics', 'Biostatistics (Michigan)'};
save('CommunityResults_firstlayer.mat', 'A','authorNames','labels_firstlayer','clusterNames_firstlayer');
                      




%%%% Generate the scree plot and the permuted adjacency matrix plot%%%%%%%%

%%% scree plot %%%%%%
[~, eigvals] = eigs(A+eye(n), 50);
eigvals = diag(eigvals);
[eigvals,~] = sort(abs(eigvals),'descend');
eigvals = abs(eigvals);
plot(eigvals,'-rs','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'scree.pdf');


%%% permuted adjacency matrix %%%%%%
[~, ix] = sort(labels_firstlayer, 'ascend');   
A_draw = A(ix,ix);    
spy(A_draw);
xlabel('');
rgb = [50 215 50]/255;
set(get(gca,'children'),'markeredgecolor',rgb);

%%% add cluster labels %%%%%%%%
summary = tabulate(labels_firstlayer);
class_size = summary(:,2);
vec = cumsum(class_size);
vec = [1; vec];
midpoint = (vec(1:end-1) + vec(2:end))/2;
hold on
for j=1:length(midpoint)
    h = rectangle('Position', ...
        [vec(j) vec(j) vec(j+1)-vec(j)*[1 1]]);
    set(h, 'EdgeColor', [1 1 1]*.3);
    set(h, 'Linewidth', 2);
    h = text(midpoint(j), midpoint(j), num2str(j));
    set(h, 'color', [0 0 0]);
    set(h, 'HorizontalAlignment', 'center');
    set(h, 'FontSize', 30*(1 + max(0, (class_size(j)./(n/K)-1)/3)));
end
hold off

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'net-adjacency.pdf');





