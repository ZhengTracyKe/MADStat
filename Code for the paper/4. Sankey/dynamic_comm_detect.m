load('output/CoauSankeyFinal.mat');

names = authorNames(V);


%%%%% Apply SCORE to the network of each period%%%%%%%%%%%%%%%%%%%%%
K_all = [3 4 3];
labels_dynamic = zeros(length(V), 3);
degrees_dynamic = zeros(length(V),3);

for t = 1:3
    K = K_all(t);
    G = CoauAdjSankey{t,2};  %%% get the giant component 
    A = CoauAdjSankey{t,1}(G, G);   %%%% get the adjacency matrix
    n = size(A,1);
    tempLabels = SCORE(A + eye(n), K);
    degs = sum(A,2);
    
    %%% record the labels into the |V|-by-3 matrix 'labels_dynamic'%%%
    [~, ind_G,ind_V] = intersect(G, V);
    labels_dynamic(ind_V, t) = tempLabels(ind_G);
    degrees_dynamic(ind_V, t) = degs(ind_G); 
    
    
    %%% plot the permuted adjacency matrix %%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1, 3, t);
    [~, ix] = sort(tempLabels, 'ascend');   
    A_draw = A(ix,ix);    
    spy(A_draw);
    xlabel('');
    rgb = [50 215 50]/255;
    set(get(gca,'children'),'markeredgecolor',rgb);
    summary = tabulate(tempLabels);
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


%%%%% give names to clusters %%%%%%%%
clunames_to_assign = {{'Out(period1)'; 'Sp(period1)'; 'Np(period1)'; 'Bio(period1)'};
                      {'Np(period2)'; 'Sp/Bay(period2)'; 'Bio(period2)'; 'Bay(period2)'};
                      {'Sp/Np(period3)'; 'Out(period3)'; 'Hd(period3)'; 'Bay(period3)'}};

clusterNames = cell([3,1]);
for t = 1:3
    clustername = clunames_to_assign{t};
    summary = tabulate(labels_dynamic(:,t));
    [~, sort_id] = sort(summary(:,2), 'descend');
    clustername(sort_id) = clustername;
    clusterNames{t} = clustername;
     %%% print the summary
    clustersize = summary(:,2);
    table(clustername, clustersize) 
end




%%%%% Output representative nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% from Period 1 to Period 2
[C,~,ic] = unique(labels_dynamic(:,1:2), 'rows');
degrees_12 = degrees_dynamic(:,1)+degrees_dynamic(:,2);
for k = 1:size(C,1)
    nodes_in_flow = find(ic==k);
    [~, sort_id] = sort(degrees_12(nodes_in_flow), 'descend'); 
    nodes_in_flow = nodes_in_flow(sort_id);
    names_in_flow = names(nodes_in_flow);
    %%% print 20 representative nodes %%%%%
    cluName_from = clusterNames{1}{C(k,1)+1};  %%% C(k,1) takes values {0,1,2,3}
    cluName_to = clusterNames{2}{C(k,2)};  %%% C(k,2) takes values {1,2,3,4}
    fprintf(['---From ', cluName_from, ' to ', cluName_to, '\n\n'])
    m = min(20, length(nodes_in_flow));
    for j = 1:m
        fprintf([names_in_flow{j}, '\n']);
    end
    fprintf('\n')
end

%%% from Period 2 to Period 3
[C,~,ic] = unique(labels_dynamic(:,2:3), 'rows');
degrees_23 = degrees_dynamic(:,2)+degrees_dynamic(:,3);
for k = 1:size(C,1)
    nodes_in_flow = find(ic==k);
    [~, sort_id] = sort(degrees_23(nodes_in_flow), 'descend'); 
    nodes_in_flow = nodes_in_flow(sort_id);
    names_in_flow = names(nodes_in_flow);
    %%% print 20 representative nodes %%%%%
    cluName_from = clusterNames{2}{C(k,1)};  %%% C(k,1) takes values {1,2,3,4}
    cluName_to = clusterNames{3}{C(k,2)+1};  %%% C(k,2) takes values {0,1,2,3}
    fprintf(['---From ', cluName_from, ' to ', cluName_to, '\n\n'])
    m = min(20, length(nodes_in_flow));
    for j = 1:m
        fprintf([names_in_flow{j}, '\n']);
    end
    fprintf('\n')
end

save('CommunityResults_sankey.mat', 'names', 'labels_dynamic', 'degrees_dynamic', 'clusterNames');