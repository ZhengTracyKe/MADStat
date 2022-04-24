load('data/YearlyCoauNetworks_4journals.mat')

%%% Construct the coauthorship network in each period 
CoauAdjSankey = cell(3,2);
startYear_all = [1975, 1995, 2005];
endYear_all = [1997, 2007, 2015];

for j = 1:3
    startY = startYear_all(j);
    endY = endYear_all(j);
    rawAdj = CoauAdjYearly_4J{startY-1974};
    for k = (startY+1):endY
        rawAdj = rawAdj + CoauAdjYearly_4J{k-1974};
    end
    A = double(rawAdj>0);
    A = triu(A,1)+tril(A, -1);
    %%% restrict to the giant component
    bins = conncomp(graph(A));
    GiantComp = find(bins==mode(bins));
    CoauAdjSankey{j,1} = A;
    CoauAdjSankey{j,2} = GiantComp;
end

%%% Find the set V = (G1 & G2)||(G2 & G3) %%%%%%%%
G1 = CoauAdjSankey{1,2};
G2 = CoauAdjSankey{2,2};
G3 = CoauAdjSankey{3,2};
V = union(intersect(G1, G2),intersect(G2, G3));

%%% Read the author names %%%%%%%%%%%%%%%%%%%%%%%%

[filename, filepath] = uigetfile('data/author_name.txt');
fullname = fullfile(filepath, filename); %%% get the full path (for fopen)
fid=fopen(fullname,'r','n','UTF-8');
data = textscan(fid,'%s','Delimiter',',');
authorNames = data{1};
fclose(fid);
clear data;

save('CoauSankeyFinal.mat', 'CoauAdjSankey', 'authorNames','V');

