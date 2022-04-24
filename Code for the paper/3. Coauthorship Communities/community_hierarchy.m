load('output/CommunityResults_firstlayer');
n = size(A,1);
K = max(labels_firstlayer);

labels_hierarchy= zeros(n, 3);
for k = 1:6
    labels_hierarchy(labels_firstlayer==k,1)=k;
end


%%% Cluster 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_sub = 3;
nodes_interest = find(labels_firstlayer==1);
[~,~,clu1labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu1labels(clu1labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [1, 3, 2];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu1labels == sort_id(k)) = mapping(k);
end
clu1labels = temp;

labels_hierarchy(labels_firstlayer==1, 2) = clu1labels(labels_firstlayer==1);

%%%% print representative nodes
fprintf('------------ Community 1 ------------ \n\n');
for k = 1:K_sub
    nodes_interest = find(clu1labels==k);
    fprintf(['Community 1-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end


%%% Cluster 11: further split %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_interest = find(clu1labels==1);
K_sub = 6;
[~,~,clu11labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu11labels(clu11labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [4, 6, 1, 2, 3, 5];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu11labels == sort_id(k)) = mapping(k);
end
clu11labels = temp;

labels_hierarchy(clu1labels==1, 3) = clu11labels(clu1labels==1);

%%%% print representative nodes
fprintf('------------ Community 1-1 ------------ \n\n');
for k = 1:K_sub
    nodes_interest = find(clu11labels==k);
    fprintf(['Community 1-1-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end




%%% Cluster 2 (no split)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_interest = find(labels_firstlayer==2);
[SgnQpvals,~,~] = SgnQ(A(nodes_interest,nodes_interest));
SgnQpvals

%%%% print representative nodes
fprintf('------------ Community 2 ------------ \n\n');
PrintNodes(A, authorNames, nodes_interest);






%%% Cluster 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_sub = 3;
nodes_interest = find(labels_firstlayer==3);
[~,~,clu3labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu3labels(clu3labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [1, 3, 2];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu3labels == sort_id(k)) = mapping(k);
end
clu3labels = temp;

labels_hierarchy(labels_firstlayer==3, 2) = clu3labels(labels_firstlayer==3);

%%%% print representative nodes
fprintf('------------ Community 3 ------------ \n\n');
for k = 1:K_sub
    nodes_interest = find(clu3labels==k);
    fprintf(['Community 3-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end





%%% Cluster 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_sub = 4;
nodes_interest = find(labels_firstlayer==4);
[~,~,clu4labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu4labels(clu4labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [1 4 3 2];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu4labels == sort_id(k)) = mapping(k);
end
clu4labels = temp;

%%%% extract the cluster not in the giant component
outliers = setdiff(nodes_interest, find(clu4labels>0));
Atemp = A(outliers, outliers);
bins = conncomp(graph(Atemp));
extracluster = outliers(bins==mode(bins));
clu4labels(extracluster) = 5;


labels_hierarchy(labels_firstlayer==4, 2) = clu4labels(labels_firstlayer==4);

%%%% print representative nodes
fprintf('------------ Community 4 ------------ \n\n');
for k = 1:(K_sub+1)
    nodes_interest = find(clu4labels==k);
    fprintf(['Community 4-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end




%%% Cluster 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_sub = 7;
nodes_interest = find(labels_firstlayer==5);
[~,~,clu5labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu5labels(clu5labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [3, 6, 2, 1, 4, 7, 5];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu5labels == sort_id(k)) = mapping(k);
end
clu5labels = temp;

labels_hierarchy(labels_firstlayer==5, 2) = clu5labels(labels_firstlayer==5);

%%%% print representative nodes
fprintf('------------ Community 5 ------------ \n\n');
for k = 1:K_sub
    nodes_interest = find(clu5labels==k);
    fprintf(['Community 5-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end




%%% Cluster 53: further split %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_interest = find(clu5labels==3);
K_sub = 2;
[~,~,clu53labels,~] = NetworkPartition(A, nodes_interest, K_sub, 1);

%%%% adjust the labeling of clusters to match the paper
summary = tabulate(clu53labels(clu53labels>0));
[~, sort_id] = sort(summary(:,2), 'descend');
mapping = [2, 1];
temp = zeros(n,1);
for k=1:K_sub
    temp(clu53labels == sort_id(k)) = mapping(k);
end
clu53labels = temp;

labels_hierarchy(clu5labels==3, 3) = clu53labels(clu5labels==3);

%%%% print representative nodes
fprintf('------------ Community 5-3 ------------ \n\n');
for k = 1:K_sub
    nodes_interest = find(clu53labels==k);
    fprintf(['Community 5-3-', num2str(k), ': ', num2str(length(nodes_interest)), ' nodes\n\n'])
    PrintNodes(A, authorNames, nodes_interest);
end







%%% Cluster 6 (no split)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_interest = find(labels_firstlayer==6);
[SgnQpvals,~,~] = SgnQ(A(nodes_interest,nodes_interest));
SgnQpvals

%%%% print representative nodes
fprintf('------------ Community 6 ------------ \n\n');
PrintNodes(A, authorNames, nodes_interest);






save('CommunityResults_hierarchy.mat', 'labels_hierarchy', 'authorNames');

%%% Print the community label of any individual of interest %%%%%%%%%%%%%%%
person = 'Bradley Efron';   
label_thisperson = labels_hierarchy(strcmp(authorNames, person),:);
s = sum(label_thisperson>0);
fprintf(['Community label of ', person, ': C'])
for l = 1:s
    fprintf(num2str(label_thisperson(l)));
    if (l<s)
        fprintf('-');
    else
        fprintf('\n');
    end
end
