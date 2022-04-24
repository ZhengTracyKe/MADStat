function [GiantComp, outliers, clusterLabels, SgnQpvals] = NetworkPartition(A, nodes_interest, K, plotting)


n = size(A,1);

%%%% restrict to the subnetwork and find the giant component 
A1 = A(nodes_interest, nodes_interest);
bins = conncomp(graph(A1));
GiantComp = nodes_interest(bins==mode(bins));  
outliers = nodes_interest(bins~=mode(bins));
A_sub = A(GiantComp, GiantComp);
n_sub = size(A_sub,1);


%%%% generate the scree plot
[~, eigVals] = eigs(A_sub + eye(n_sub), 20);
eigVals = diag(eigVals);
[~, ix] = sort(abs(eigVals),'descend');
eigVals = eigVals(ix);
ind = find(eigVals<0);
if (plotting)
    subplot(1,2,1)
    plot(abs(eigVals), '-sb', 'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','b');
    hold on
    plot(ind, abs(eigVals(ind)), 'sr', 'MarkerSize',10.1, 'MarkerFaceColor','r');
    hold off
end


%%%% community detection by SCORE
tempLabels = SCORE(A_sub + eye(n_sub), K);
summary = tabulate(tempLabels);
clusterLabels = zeros(n,1);
clusterLabels(GiantComp) = tempLabels; %%%(use the index of nodes in the origial network)
if (plotting)
    subplot(1,2,2)
    [~, ix] = sort(tempLabels, 'ascend');   
    A_draw = A_sub(ix,ix);    
    spy(A_draw);
    xlabel('');
    rgb = [50 215 50]/255;
    set(get(gca,'children'),'markeredgecolor',rgb);
    %%% add cluster labels %%%%%%%%
    class_size = summary(:,2);
    vec = [1; cumsum(class_size)];
    midpoint = (vec(1:end-1) + vec(2:end))/2;
    hold on
    for j=1:length(midpoint)
        h = rectangle('Position', [vec(j) vec(j) vec(j+1)-vec(j)*[1 1]]);
        set(h, 'EdgeColor', [1 1 1]*.3);
        set(h, 'Linewidth', 2);
        h = text(midpoint(j), midpoint(j), num2str(j));
        set(h, 'color', [0 0 0]);
        set(h, 'HorizontalAlignment', 'center');
        set(h, 'FontSize', 30*(1 + max(0, (class_size(j)./(n/K)-1)/3)));
    end
    hold off
end


%%%% compute the SgnQ pvalues for each community
SgnQpvals = zeros(K,1);
for k = 1:K
    A_bycluster = A(clusterLabels==k,clusterLabels==k);
    [SgnQpvals(k),~,~] = SgnQ(A_bycluster);
end


%%%% printing the summary 
Cluster = summary(:,1);
Size = summary(:,2);
Percent = summary(:,3);
table(Cluster, Size, Percent, SgnQpvals)



end