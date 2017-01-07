function d=sqdistALL(a,b)
% SQDIST - computes squared Euclidean distance matrix
%          computes a rectangular matrix of pairwise distances
% between points in A  and points in B

% NB: very fast implementation taken from Roland Bunschoten

%modified 21 may 2014
%-----------------
a=a';
b=b';
%-----------------
 
aa = sum(a.*a,1); bb = sum(b.*b,1); ab = a'*b; 
d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);