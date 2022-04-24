load('data/YearlyCiteeNetworks.mat');   
%%%this file contains the ID of 2831 nodes in the citee network, stored
%%%in the variable "keepNodeID"

load('data/CiteeAdjFinal.mat');



%%% aggregate the yearly citee networks %%%%%
startYear_all = [1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011];
endYear_all = [2000,2001,2001,2002,2003,2004,2004,2005,2006,2007,2007,2008,2009,2010,2010,2011,2012,2013,2013,2014,2015];
k0 = 1974;

CiteeAdjAggregate = cell(21,1);

for j = 1:21
    startY = startYear_all(j);
    endY = endYear_all(j);
    rawAdj = CiteeAdjYearly{startY-1974};
    for k = (startY+1):endY
        rawAdj = rawAdj + CiteeAdjYearly{k-1974};
    end
    A = double(rawAdj(keepNodeID,keepNodeID)>=2);
    A = triu(A,1)+tril(A, -1);
    CiteeAdjAggregate{j} = A;
end
%%% Note: Each aggregrated network is restricted to the 2831 nodes


save('CiteeDynamicFinal.mat','CiteeAdjAggregate','keepNodeID');