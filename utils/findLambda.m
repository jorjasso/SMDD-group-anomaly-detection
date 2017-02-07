function [Q1,Q3,Q5,Q7,Q9]=findLambda(S)
% function lambda=findLambda(S); finds the lambda parameter for a kernel of
% probability measures given by the RBF kernel based on the quantile
% function

if size(S,1)<size(S,2)
    S=S';
end
% huge dataset = number of samples*number of dimensions
factor=sum(cell2mat(cellfun(@(x) size(x,1),S,'UniformOutput',0)))*size(S{1},2)

if factor <10000
    
    tic;
    X=cell2mat(S);
    M=sqdistAll(X,X);
    Q1=quantile(M(:),0.1)
    Q3=quantile(M(:),0.3)
    Q5=quantile(M(:),0.5)
    Q7=quantile(M(:),0.7)
    Q9=quantile(M(:),0.9)
    median(M(:))
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
    M=sqdistAll(X,X);
    Q1=quantile(M(:),0.1)
    Q3=quantile(M(:),0.3)
    Q5=quantile(M(:),0.5)
    Q7=quantile(M(:),0.7)
    Q9=quantile(M(:),0.9)
end

