function label = SCORE(A, K, truncate, t)

%%%INPUT%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% A:         adjacency matrix of a symmetric network
%%%
%%% K:         number of communities
%%%
%%% truncate:  if this is a nonzero number, truncate the eigen-ratio matrix;
%%%            default value is 1;
%%%            set it to zero if truncation is not wanted
%%%
%%% t:         the truncation threshold, when 'truncate' is a nonzero number;
%%%            default value is log(n)
%%%
%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% label:     the community lables of all nodes
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(A,1);	                % Number of nodes

% Initical checks
if (~issymmetric(A))
    fprint('Error: The adjacency matrix is not symmetrical.\n Please use the version of SCORE for directed networks.')
    return
end
if (max(conncomp(graph(A)))>1)
    fprintf('Error: The network is not connected.\n Please apply SCORE to the giant component.\n');
    return  
end
if K < 2 
    fprintf('Error: K should be larger than 1.\n');  
    return
end

% Set the threshold values
if (nargin ==2)
    truncate = 1;
    t = log(n);
end
if (nargin ==3)
    if (truncate~=0)
        t = log(n);
    end
end


% SCORE
[vb,~] = eigs(A,K);
vb2 = vb(:,2:K);
rb = vb2./repmat(vb(:,1),[1,K-1]);
if (truncate)
    rb(rb > t) = t;   
    rb(rb < -t) = -t;
end
label = kmeans(rb,K,'Replicates',1000, 'MaxIter', 200);
end
