function SelectedNodes = PrintNodes(A,names, nodes_interest, howMany)


%%%% By default, output 20 nodes with highest degrees, 5 nodes with highest
%%%% betweenness and 5 nodes with highest closedness 
if (nargin<4)
    howMany = [20, 5, 5]; 
end
N_deg = howMany(1);
N_btw = howMany(2);
N_cls = howMany(3);



A_sub = A(nodes_interest, nodes_interest);
G = graph(A_sub);

degrees = centrality(G, 'degree');
[~, sort_id] = sort(degrees, 'descend');
HighDegNodes = nodes_interest(sort_id(1:N_deg)); 
fprintf('High degree nodes:\n');
for k = 1:N_deg
    fprintf([names{HighDegNodes(k)},',\n']);
end
fprintf('\n');


betweenness = centrality(G, 'betweenness');
[~, sort_id] = sort(betweenness, 'descend');
HighBtwNodes = nodes_interest(sort_id(1:N_btw)); 
fprintf('High betweenness nodes: ');
for k = 1:N_btw
    fprintf([names{HighBtwNodes(k)},', ']);
end
fprintf('\n');

closeness = centrality(G, 'closeness');
[~, sort_id] = sort(closeness, 'descend');
HighClsNodes = nodes_interest(sort_id(1:N_cls)); 
fprintf('High closeness nodes: ');
for k = 1:N_cls
    fprintf([names{HighClsNodes(k)},', ']);
end
fprintf('\n\n');

SelectedNodes = {HighDegNodes,HighBtwNodes,HighClsNodes};

end