function [lambda firstQuantile lastQuantile]=findLambda(S)
% function lambda=findLambda(S); finds the lambda parameter for a kernel of
% probability measures given by the RBF kernel
%
% tic
% n=length(S);
% val=cell(n,n);
% sum=0;
% for i=1:length(S)
%     for j=1:length(S)
%         X=S{i}; Z=S{j};
%         dX=diag(X*X'); dZ=diag(Z*Z');
%         XX=dX*ones(1,size(Z,1));ZZ=ones(size(X,1),1)*dZ';
%         val{i,j}=XX-2*X*Z'+ZZ ;
%         sum=sum+length(val{i,j}(:));
%     end
% end
%
% M=cell2mat(val);
% lambda=median(M(:));
% firstQuantile=quantile(M(:),0.1);
% lastQuantile=quantile(M(:),0.9);
% toc
% [lambda firstQuantile lastQuantile]
%another code
%standarize input
if size(S,1)<size(S,2)
    S=S';
end
% big data?= number of samples*number of dimensions
factor=sum(cell2mat(cellfun(@(x) size(x,1),S,'UniformOutput',0)))*size(S{1},2)

if factor <100000
    
    tic;
    X=cell2mat(S);
  %  Z=X;
  %  dX=diag(X*X'); dZ=diag(Z*Z');
  %  XX=dX*ones(1,size(Z,1));ZZ=ones(size(X,1),1)*dZ';
  %  M=XX-2*X*Z'+ZZ;
     M=sqdistAll(X,X);
    lambda=median(M(:));
    firstQuantile=quantile(M(:),0.1);
    lastQuantile=quantile(M(:),0.9);
    toc        
else
    % take randomly 10000 points
    
    %aproximate code
    disp('aproximate lambda computation')
    tic
    nroSamples=10;
    %minimum samples accross matrices in the cell
    val=min(cell2mat(cellfun(@(x) size(x,1),S,'UniformOutput',0)));
    if val<10
        nroSamples=val;
    end
    T=cell(length(S),1);
    for i=1:length(S)
        % ten (nroSamples) randomly choosen by group
        indImage=randperm(size(S{i},1),nroSamples);
        T{i}=S{i}(indImage,:);
    end
    
    X=cell2mat(T);
    %Z=X;
    %dX=diag(X*X'); dZ=diag(Z*Z');
    %XX=dX*ones(1,size(Z,1));ZZ=ones(size(X,1),1)*dZ';
    %M=XX-2*X*Z'+ZZ;
    M=sqdistAll(X,X);
    lambda=median(M(:));
    firstQuantile=quantile(M(:),0.1);
    lastQuantile=quantile(M(:),0.9);
    toc    
        
end

