function experimentGroupAnomalyDetectionGMM(option,nroEpochs)
%   experimentGroupAnomalyDetectionGMM(option,nroEpochs) 
%         option ={0,2} Distribution-based group anomaly detection
%         experiments   over a GMM dataset
%         option=1 Point-based group anomaly detection  over a GMM dataset, 
%         nroEpochs= number of runs

run   /home/jorjasso/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
addpath ./SVM-KM/


nTraining=50;
nTest=30;
kernelOp=3;
kappa=ones(nTraining,1);
sigma1Values=zeros(1,200);sigma2Values=zeros(1,200);sigma3Values=zeros(1,200);
yt=[ones(20,1);-ones(10,1)]; %there are 20 anomalous groups and 10 normal groups for test
for epoch=1:nroEpochs
    
    switch option
        case 0 % distribution based group anomalies
             [training,  test]=getDistributionBasedData(nTraining,nTest,0);
        case 1 % point based group anomalies
            [training,  test]=getPointBasedData(nTraining,nTest,0);
        case 2 % distribution based group anomalies, another case
            [training,  test]=generateAnomalusGroupDetectionData(nTraining,nTest,0);
        otherwise
            disp('no valid option')
    end
    if option==0
       
        
    else
        
    end
    
    %kernel parameter
    S=training{1};
    
    [sigma1 sigma2 sigma3 ]=findLambda(S);
    
    sigma1Values(epoch)=sigma1;
    sigma2Values(epoch)=sigma2;
    sigma3Values(epoch)=sigma3;
    
    % inputs
    S=training{1}; mu=training{2};
    STest=test{1}; muTest=test{2};
    
    %Experiment
    %----------------
    k=1;
    for kernelParam=[1/sigma1 1/sigma2 1/sigma3]
        % KERNEL MATRICES
        %------------------------
        % kernels for SMDDCPP, SMDDDA and OCSMM
        %-------------------------------------------
        disp('computing kernels')
        tic
        [Krr,Kre,Kee,K_traceTraining,K_traceTest]=kernelPM(S, STest,kernelOp,kernelParam); 
        
        % kernels for SMDDA with spherical normalization
        %-------------------------------------------
         [Krr_SN,Kre_SN,Kee_SN,Borrar]=kernelPM(S, STest,4,kernelParam);
         
        %kernels for SVDD
        %-------------------------------------------
        Krr_SVDD=computeKernel(mu,mu,kernelOp,kernelParam);
        Kre_SVDD=computeKernel(mu,muTest,kernelOp,kernelParam);
        if kernelOp==3% RBF kernel
            Kee_SVDD=ones(size(muTest,1),1);
        else
            Kee_SVDD=diag(computeKernel(muTest,muTest,kernelOp,kernelParam))';
        end
        toc
        %------------------------
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
        k=k+1;
    end
    %----------------------------------------------------------------------   

end

switch option
    case 0
        str='Unimodal'
    case 1
        str='Multimodal'
    case 2
        str='Another'
end

save (['experimentosGroupAnomaluGMM' str '.mat'])
%load (['experimentosGroupAnomaluGMM' str '.mat'])

% choose regularization and RBF paramters:
% for example (C=1,lamda=2) correspond (C=1, 1/sigma2)

C=1;lambda=1
matrixC = cell2mat(statistics(:,lambda,:,C)); %[epoch ModelOp C]
%matrixC = cell2mat(statistics(:,:,C,lambda)); %[epoch ModelOp C]

%plots statistics
ind=mod(1:50,10);
figure
boxplot(100-matrixC(:,ind==1),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline');           %ACC
ylabel('ACC')
title ('Accuracy')
figure
boxplot(matrixC(:,ind==2),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline');               %AUC
ylabel('AUC')
title ('Area Under the Curve Roc')
figure
boxplot(100-matrixC(:,ind==3),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % ACC anomalous groups
ylabel('ACC ')
title ('Accuracy Anomalous groups')
figure
boxplot(100-matrixC(:,ind==4),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % ACC normal groups
ylabel('ACC ')
title ('Accuracy Normal groups')
figure
boxplot(matrixC(:,ind==5),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % number of support vector expansion
ylabel('nro SV ')
title ('Number of support vectors')
% figure
% boxplot(matrixC(:,ind==6),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % number of support vectors
figure
boxplot(matrixC(:,ind==7),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % number of support vectors errors
ylabel('nro SV ')
title ('Number of support vectors errors')
figure
boxplot(matrixC(:,ind==8),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % C
ylabel('[0 1]')
title ('Regularization Parameter C')
figure
boxplot(matrixC(:,ind==9),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % F1
ylabel('F1 ')
title ('F1 score')
figure
boxplot(matrixC(:,ind==0),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline'); % MCC
ylabel('MCC')
title ('Matthews correlation coefficient')
figure
boxplot(1./eval(['sigma' int2str(lambda) 'Values'])); %lambda
ylabel('lambda')
title ('RBF kernel parameter lambda')

return


%PRINTS
%---------------------------
%AUC
set(gca,'FontSize',30,'fontWeight','bold')
set(gca,'XTickLabel',{' '})
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
%boxplot(matrixAUC(:,ind==2),'labels',{'SMDDCPP','SMDDDA','SMDDDASN','OCSMM','SVDD'},'labelorientation','inline');
boxplot(matrixC(:,ind==2));
print(gcf, '-dpng', '-r0', ['GMM/AUC2.png']);

%ACC anomalous
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
boxplot(100-matrixC(:,ind==3)); %lambda
print(gcf, '-dpng', '-r0', ['GMM/ACCA2.png']);

%ACC normal
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
boxplot(100-matrixC(:,ind==4)); %lambda
print(gcf, '-dpng', '-r0', ['GMM/ACCN2.png']);

%Print  groups
%-------------
if option==0
    [training,  test]=getDistributionBasedData(nTraining,nTest,0);
else
    [training,  test]=getPointBasedData(nTraining,nTest,0);
end

%print means of groups
S=training{1};muTrain=training{2};Sigma=training{3};
SVal=test{1};muVal=test{2};SigmaTest=test{3}


figure
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
plot(muTrain(:,1),muTrain(:,2),'.g','MarkerSize',30) % means training set
hold on
plot(muVal(1:20,1),muVal(1:20,2),'.r','MarkerSize',30)  % means group of anomalous points
%plot(muVal(11:20,1),muVal(11:20,2),'.r','MarkerSize',30) % means of group anomalies
%plot(muVal(21:30,1),muVal(21:30,2),'.k','MarkerSize',30) % means of normal
%groups
saveas(gcf,['GMM/meanDistribution'],'fig');
print(gcf, '-dpng', '-r0', ['GMM/meanDistribution.png']);


%for sloan
%AUC
set(gca,'FontSize',30,'fontWeight','bold')
%set(gca,'XTickLabel',{' '})
boxplot(matrixC(:,ind==2),'labels',{'M1','M2','M3','M4','M5'});
%set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
set(findobj(gca,'Type','text'),'FontSize',30)
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
print(gcf, '-dpng', '-r0', ['sloan/AUC2_4.png']);

%ACC anomalous
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
boxplot(100-matrixC(:,ind==3),'labels',{'M1','M2','M3','M4','M5'}); %lambda
set(findobj(gca,'Type','text'),'FontSize',30)
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
print(gcf, '-dpng', '-r0', ['sloan/ACCA2_4.png']);

%ACC normal
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
boxplot(100-matrixC(:,ind==4),'labels',{'M1','M2','M3','M4','M5'}); %lambda
set(findobj(gca,'Type','text'),'FontSize',30)
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
print(gcf, '-dpng', '-r0', ['sloan/ACCN2_4.png']);

tam=30
figure; % solo las medias
subplot(2,2,1)
set(gca,'FontSize',tam,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',tam,'fontWeight','bold')
%plot(mu(:,1),mu(:,2),'.g');
scatter(mu(:,1),mu(:,2),300,'g.')
hold on;
%plot(muTest(1:50,1),muTest(1:50,2),'.r');
scatter(muTest(1:50,1),muTest(1:50,2),300,'.r')
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
subplot(2,2,2)
set(gca,'FontSize',tam,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',tam,'fontWeight','bold')
scatter(mu(:,2),mu(:,3),300,'g.')
hold on;
%plot(muTest(1:50,1),muTest(1:50,2),'.r');
scatter(muTest(1:50,2),muTest(1:50,3),300,'.r')
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
subplot(2,2,3)
set(gca,'FontSize',tam,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',tam,'fontWeight','bold')
scatter(mu(:,3),mu(:,4),300,'g.')
hold on;
%plot(muTest(1:50,1),muTest(1:50,2),'.r');
scatter(muTest(1:50,3),muTest(1:50,4),300,'.r')
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
subplot(2,2,4)
set(gca,'FontSize',tam,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',tam,'fontWeight','bold')
scatter(mu(:,4),mu(:,1),300,'g.')
hold on;
%plot(muTest(1:50,1),muTest(1:50,2),'.r');
scatter(muTest(1:50,4),muTest(1:50,1),300,'.r')
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
print(gcf, '-dpng', '-r0', ['sloan/meanDistribution_4.png']);
