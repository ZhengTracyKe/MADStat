function Adj = FuncGetNetwork(AuPap, PapPap, type, startYear, endYear, citLength, selfCite)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AuPap: the Author-Paper list (columns: authorID, paperID, paperYear,
% paperJournal
%
% PapPap: the Paper-Paper-Citaion list (columns: citePaperID, citedPaperID,
% citePaperYear, citedPaperYear, Is_selfcite)
%
% Type: the typpe of networks to construct (values: citee, citer, citation)
%
% [startYear, endYear]: the range of citePaperYear
%
% citLength: the maximum difference between citePaperYear and
% citedPaperYear
% 
% selfCite: 0 means excluding self citations, 1 means including self
% citations
%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nau = max(AuPap(:,1));
Npap = max(AuPap(:,2));

if nargin < 7
    selfCite = 0;   % default: not including self citations
else if nargin < 6
        citLength = Inf;   % defacult: no restriction on citation length
    else if nargin < 5
            endYear = max(AuPap(:,3));   %default: no restriction on end year
        else if nargin < 4
            startYear = max(AuPap(:,3));   %default: no restriction on start year
            end
        end
    end
end

%%%% determine the reference records to retain %%%%%%%%%%%%%%%%%%%%%%%
keepInd = (PapPap(:,3)-PapPap(:,4)<citLength & PapPap(:,3)-PapPap(:,4)>=0 ...
    & PapPap(:,3)>=startYear & PapPap(:,3)<=endYear);
if selfCite==0
    keepInd = (keepInd & PapPap(:,5)==0);
end
PapPap_new = PapPap(keepInd,:);


%%%% construct the adjacency matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = sparse(AuPap(:,1),AuPap(:,2),1,Nau,Npap);
Y = sparse(PapPap_new(:,1),PapPap_new(:,2),1,Npap,Npap);

if strcmp(type,'citee')
    Adj = X*(Y'*Y)*X';
end

if strcmp(type,'citer')
    Adj = X*(Y*Y')*X';
end

if strcmp(type, 'citation')
    Adj = X*Y*X';
end




end