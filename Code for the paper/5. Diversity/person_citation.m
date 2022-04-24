clear;
load('data/AuthorPaperMatrices.mat')

%%
Nauthor = max(AuPap(:,1));
Npaper = max(AuPap(:,2));
Year_all = 1975:2015;

Adj = sparse(Nauthor,Nauthor);

for k = 1:length(Year_all)
    Adj = Adj + FuncGetNetwork(AuPap, PapPap, 'citation', ...
        Year_all(k), Year_all(k), 10, 0);
end

%% node screening
A = Adj;
A = double(A > 0);
A = triu(A,1) + tril(A,-1);

clear AuPap;
clear PapPap;
clear X;
clear Y;
clear Adj;
clear k;

[filename, filepath] = uigetfile('data/author_name.txt');
fullname = fullfile(filepath, filename); %%% get the full path (for fopen)
fid=fopen(fullname,'r','n','UTF-8');
data = textscan(fid,'%s','Delimiter',',');
name = data{1};
fclose(fid);
clear data;
clear filename;
clear filepath;
clear fullname;
clear fid;

%% high-degree nodes
deg = sum(A);
[~, ranking] = sort(deg,'descend');                 % get ranking by degree
t1 = 1;
t2 = 20;
top_name = name(ranking(t1:t2));     
top_name

%% the 1000-th biggest degree
deg(ranking(1000))
%% top 100 pvs
clc;
top = 1000;
Q_pv_list_citee = zeros(top,1);
Q_st_list_citee = zeros(top,1);
N_net_list_citee = zeros(top,1);

Q_pv_list_citer = zeros(top,1);
Q_st_list_citer = zeros(top,1);
N_net_list_citer = zeros(top,1);

echo = 1;
for t = 1:top
    name_hub = name{ranking(t)}; 
    id_hub = find(strcmp(name,name_hub));
    id_hub = id_hub(1);
    
    id_citee = find(A(id_hub,:)~=0);
    A_citee = A(id_citee,id_citee);
    A_citee = A_citee' + A_citee;
    A_citee = double(A_citee > 0);
    
    N_net_list_citee(t) = size(A_citee,1);
    [Q_pv_list_citee(t),tmp,Q_st_list_citee(t)] = SgnQ(A_citee);
    
    id_citer = find(A(:,id_hub)~=0);
    A_citer = A(id_citer,id_citer);
    A_citer = A_citer' + A_citer;
    A_citer = boolean(A_citer > 0) * 1;
    
    N_net_list_citer(t) = size(A_citer,1);
    [Q_pv_list_citer(t),tmp,Q_st_list_citer(t)] = SgnQ(A_citer);
    
    if echo == 1
        fprintf([pad(name_hub,28),'\n']);
        fprintf(pad(name_hub,28));
        fprintf([' | network size: ',pad(num2str(N_net_list_citee(t)),10)]);
        fprintf([' | SgnQ p-value: ',pad(num2str(Q_pv_list_citee(t)),10)]);
        fprintf([' | Q-stat-citee: ',pad(num2str(Q_st_list_citee(t)),10),'\n']);
%         fprintf([' | Q-stat-citer: ',pad(num2str(Q_st_list_citer(t)),10),'\n']);
    end
end
%%
t = 2;
Q_st_list_citee(t)
Q_st_list_citer(t)
%%
name_ranked = name(ranking);
save('output/Q_stat_citee_citer.mat','name','A','name_ranked','ranking','N_net_list_citee','N_net_list_citer','Q_pv_list_citee','Q_pv_list_citer','Q_st_list_citee','Q_st_list_citer');