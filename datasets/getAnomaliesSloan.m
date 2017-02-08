function [nonAnomalous, anomalousFirstType,anomalousSecondType,anomalousThirdType] = getAnomaliesSloan(N)
%  This function get the non anomalous galaxy clusters and N non anomalous groups, (three different settings)
%

%%%%%%%%%%%%%%%%%
%NON-ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%

%reading data
%-------------
M = csvread('spectrum_s1.csv',1,0);
ind=mod(1:15059,2);
MM=M(ind==1,:);
ind=MM(:,1); %indices of clusters

% scaling data to [-1 1]
%-------------
X=bsxfun(@rdivide, 2*bsxfun(@minus, MM(:,2:end),min(MM(:,2:end))), (max(MM(:,2:end))-min(MM(:,2:end))))-1;

%PCA
%-------------
[~,score,latent,~] = princomp(X);
AAA=latent*100/sum(latent);
sum(AAA(1:4)) % four dimensions preserves about 85% of the variance
XPCA=score(:,1:4);
%free memory
%-------------
clear M MM X;

%clusters of galaxies as cells, computing the mean of the covariance of all
%clusters
%-------------------------------------------------------------------------
temp=0;

S=[];mu=[];Sigma=[];
for i=0:504
    S{i+1}=XPCA(ind==i,:);
    temp=temp+cov(S{i+1});
    mu(i+1,:)=mean(S{i+1});Sigma{i+1}=cov(S{i+1});
end
SigmaMean=temp/505; %expectation of the covariance matrix  from  all nonanomalous groups

nonAnomalous={S mu Sigma};

%%%%%%%%%%%%%%%%%
%ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%
%First type of group anomalies
S=[];mu=[];Sigma=[];
[m,~]=size(XPCA);

[minVal,~]=size(XPCA(ind==1,:));
[maxVal,~]=size(XPCA(ind==504,:));
    
for i=1:N
  
    nroGalaxies=randi([maxVal,minVal],1,1);
    n=poissrnd(nroGalaxies);
    r = randi(m,n,1);
    S{i}=XPCA(r,:);
    mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
anomalousFirstType={S mu Sigma};

%Second type of group anomalies
S=[];mu=[];Sigma=[];
[m,~]=size(XPCA);
for i=1:N
    
    values=randi([maxVal,minVal],3,1);
   
    % group anomalies generated from a GMM
    n=[poissrnd(values(1));poissrnd(values(2));poissrnd(values(3))];
    
    muGM=[mean(XPCA(randi(m,n(1),1),:)); mean(XPCA(randi(m,n(2 ),1),:)); mean(XPCA(randi(m,n(3),1),:))];
    sigmaGM = 0.5*SigmaMean;
    p=[0.33,0.33,0.33];%weigths for normal groups
    obj = gmdistribution(muGM,sigmaGM,p);
    nroGalaxies=randi([maxVal,minVal],1,1);
    S{i} = random(obj,poissrnd(nroGalaxies));
    mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
anomalousSecondType={S mu Sigma};

%Third type of group anomalies

for i=1:N
    
    values=randi([maxVal,minVal],3,1);
   
    % group anomalies generated from a GMM
    n=[poissrnd(values(1));poissrnd(values(2));poissrnd(values(3))];
    
    muGM=[mean(XPCA(randi(m,n(1),1),:)); mean(XPCA(randi(m,n(2 ),1),:)); mean(XPCA(randi(m,n(3),1),:))];
    sigmaGM = SigmaMean;
    p=[0.33,0.33,0.33];%weigths for normal groups
    obj = gmdistribution(muGM,sigmaGM,p);
    nroGalaxies=randi([maxVal,minVal],1,1);
    S{i} = random(obj,poissrnd(nroGalaxies));
    mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
anomalousThirdType={S mu Sigma};

end

