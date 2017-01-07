function experimentSloan(factor,nRuns)
% Experiments for the Sloan sky survey dataset (SDSS)
% factor*SigmaMean= used to construct anomalous groups from a GMM model
% with mean, and covariance = SigmaMean, SigmaMean is the mean of the
% covariance matrix of all the groups.
% If factor = 0, then each anomalous group is the cluster of randomly
% chossen galaxies if the dataset


%read data
%run   /home/jorjasso/cvx/cvx_startup.m
%addpath ./SVM-KM/
%run /home/jorjasso/Downloads/cvx/cvx_startup.m

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

%free memory
%-------------
clear M MM;

%ploting
%-------------
figure
plot(X(:,1),X(:,2),'.r');
figure
XPCA=score(:,1:4);
plot(XPCA(:,1),XPCA(:,2),'.r');

save output/sloanCompletePca.mat

load output/sloanCompletePca
%free memory
%-------------
clear X;

%clusters of galaxies as cells, computing the mean of the covariance of all
%clusters
%-------------------------------------------------------------------------
temp=0;
SNormal=cell(1,505);
for i=0:504
    SNormal{i+1}=XPCA(ind==i,:);
    temp=temp+cov(SNormal{i+1});
end
SigmaMean=temp/505; %mean cov of all normal groups


%options
[m,~]=size(XPCA);

nTraining=455;

kernelOp=3;
kappa=ones(nTraining,1);
sigma1Values=zeros(1,200);sigma2Values=zeros(1,200);sigma3Values=zeros(1,200);
yt=[ones(50,1);-ones(50,1)]; %there are 50 anomalous groups and 50 normal groups for test

%200 runs to get statistics
%------------
for epoch=1:nRuns % t
    
    %randomly form 50 anomalous groups
    SAnomalous=cell(1,50);
    for i=1:50
        if factor==0
            % unifrom random selection of galaxies to form anomalous
            % clusters
            n=poissrnd(15);            %number of points per group
            r = randi(m,n,1);
            SAnomalous{i}=XPCA(r,:);
        else
            % group anomalies generated from a GMM
            n=[poissrnd(15);poissrnd(15);poissrnd(15)];
            muGM=[mean(XPCA(randi(m,n(1),1),:)); mean(XPCA(randi(m,n(2 ),1),:)); mean(XPCA(randi(m,n(3),1),:))];
            %three anomalous topics
            %I=0.2*eye(4);
            sigmaGM = factor*SigmaMean;
            p=[0.33,0.33,0.33];%weigths for normal groups
            obj = gmdistribution(muGM,sigmaGM,p);
            SAnomalous{i} = random(obj,poissrnd(15));
        end
        
    end
    
    %free memory
    %-------------
    clear n r muGM sigmaGM p obj;
    
    %randomly separate the training  and test set
    % 455 clusters about 82% of total of clusters :82/100*(555)
    indTr = randsample(505,nTraining);
    % 50 clusters about 10% of normal clusters:  10/100*(505)
    
    indTest=setdiff([1:1:505],indTr);
    
    %training set: about 82% of total of clusters :82/100*(555)
    S=SNormal(indTr);
    %test set: about 18%  of total of clusters :82/100*(555)
    STest=[SAnomalous  SNormal(indTest)];
    
    %kernel parameters: median, 0.1 and 0.9 quantile of the distances on
    %data
    [sigma1 sigma2 sigma3 ]=findLambda(S);
    sigma1Values(epoch)=sigma1;
    sigma2Values(epoch)=sigma2;
    sigma3Values(epoch)=sigma3;
    
    
    muCell=cellfun(@mean ,S,'UniformOutput',0);
    CC=cellfun(@transpose ,muCell,'UniformOutput',0);
    mu=cell2mat(CC)';
    Sigma=cellfun(@cov ,S,'UniformOutput',0);
    
    muCellT=cellfun(@mean ,STest,'UniformOutput',0);
    CCT=cellfun(@transpose ,muCellT,'UniformOutput',0);
    muTest=cell2mat(CCT)';
    SigmaTest=cellfun(@cov ,STest,'UniformOutput',0);
    
    training={S mu Sigma};
    test={STest muTest SigmaTest};
    
    %     figure; % solo las medias
    %     %hold on; plot(score(:,1),score(:,2),'.r');
    %     hold on; plot(mu(:,1),mu(:,2),'.r');
    %     hold on;plot(muTest(51:100,1),muTest(51:100,2),'.g');
    %     hold on;plot(muTest(1:50,1),muTest(1:50,2),'.k');
    %
    %     figure; % todos los datos, vs las medias de los datos de test
    %     hold on; plot(score(:,1),score(:,2),'.r');
    %     %hold on; plot(mu(:,1),mu(:,2),'.r');
    %     hold on;plot(muTest(51:100,1),muTest(51:100,2),'.g');
    %     hold on;plot(muTest(1:50,1),muTest(1:50,2),'.k');
    
    %Experiment
    %----------------
    k=1;
    for kernelParam=[1/sigma1 1/sigma2 1/sigma3]
        % KERNEL MATRICES
        %------------------------
        % kernels for SMDDCPP, SMDDDA and OCSMM
        %-------------------------------------------
        disp('computing kernels')
        disp('kernel PM')
        tic
        [Krr,Kre,Kee,K_traceTraining,K_traceTest]=kernelPM(S, STest,kernelOp,kernelParam);
        toc
        % kernels for SMDDA with spherical normalization
        %-------------------------------------------
        disp('spherical normalization')
        tic
        [Krr_SN,Kre_SN,Kee_SN,~]=kernelPM(S, STest,4,kernelParam);
        toc
        
        %kernels for SVDD
        %-------------------------------------------
        disp('SVDD kernels')
        tic
        Krr_SVDD=kernel(mu,mu,kernelOp,kernelParam);
        Kre_SVDD=kernel(mu,muTest,kernelOp,kernelParam);
        if kernelOp==3% RBF kernel
            Kee_SVDD=ones(size(muTest,1),1);
        else
            Kee_SVDD=diag(kernel(muTest,muTest,kernelOp,kernelParam))';
        end
        toc
        %------------------------
        disp('experiments in all Models')
        tic
        for ModelOp=[1:5] %for SMDDCCP, SMDDDA SMDDDASN OCSMM SVDD
            j=1;
            for C=[1 1./(nTraining*[0.1:0.1:0.2])]
                [epoch k ModelOp j ]
                
                if (ModelOp==1||ModelOp==2||ModelOp==4)
                    kernelMatrices={Krr,Kre,Kee,K_traceTraining,K_traceTest};
                    training=S; test=STest;
                end
                if ModelOp==3
                    kernelMatrices={Krr_SN,Kre_SN,Kee_SN};
                    training=S; test=STest;
                end
                if ModelOp==5
                    kernelMatrices={Krr_SVDD,Kre_SVDD,Kee_SVDD};
                    training=mu; test=muTest;
                end
                statistics{epoch,k,ModelOp,j}=prediction_MONQP(kernelMatrices, training, test,C,kappa,kernelOp,ModelOp,yt);
                j=j+1;
            end
        end
        toc
        k=k+1;
    end
   
end

save (['output/workspaceVariablesExperimentosSloan' int2str(factor) 'Sigma.mat'])
% parece que va a demorar 11 y 13 dias el experimento