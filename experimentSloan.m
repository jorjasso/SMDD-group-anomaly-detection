function experimentSloan(option)
% Experiments for the Sloan sky survey dataset (SDSS)
%run /Users/jorgeluisguevaradiaz/Documents/GITProjects/cvx/cvx_startup.m
run   /home/jorjasso/cvx/cvx_startup.m
addpath ./SVM-KM/
addpath ./datasets/
addpath ./experiments/
addpath ./models/
addpath ./kernels/
addpath ./utils/

%%%%%%%%%%%%%%%%%%%%%%%%
%EXPERIMENT DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%
N=505
kFold=5;
kernelOp=3
%%%%%%%%%%%%%%%%%%%%%%%%
%DATA
%%%%%%%%%%%%%%%%%%%%%%%%
switch option
    case 1 % point-based group anomalies first type
        type='SloanRandomFirstTypeAnomalies';
        [nonAnomalous, anomalous,~,~] = getAnomaliesSloan(N);
    case 2 % point-based group anomalies second type
        type='SloanSecondTypeAnomalies';
        [nonAnomalous, ~,anomalous,~] = getAnomaliesSloan(N);
    case 3 % point-based group anomalies third type
        type='SloanThirdTypeAnomalies';
        [nonAnomalous, ~,~,anomalous] = getAnomaliesSloan(N);
    otherwise
        disp('no valid option')
end

%percent anomalies
nAnomalies=50;
nNonAnomalies=N;

%clases
y=[ones(nNonAnomalies,1);-ones(nAnomalies,1)];

%kappa values for SMDD CCP
kappa=ones(length(y),1);


%index: selection nNonanomalies groups from nonAnomalous dataset and nAnomalies
%groups from anomalous dataset
index=randperm(N,nNonAnomalies);
NA={nonAnomalous{1}(index),nonAnomalous{2}(index,:),nonAnomalous{3}(index)};
index=randperm(N,nAnomalies);
A={anomalous{1}(index),anomalous{2}(index,:),anomalous{3}(index)};

%Combined dataset
tempData=[NA; A ];

dataset{1}=[tempData{1,1} tempData{2,1}];  %groups
dataset{2}=[tempData{1,2} ;tempData{2,2}]; %means
dataset{3}=[tempData{1,3} tempData{2,3}];  %coovariance matrices
dataset{4}=y;                           %classes
dataset{5}=kappa;                         %kappa values

%%%%%%%%%%%%%%%%%%%%%%%%
%Nested cross validation
%%%%%%%%%%%%%%%%%%%%%%%%
CVP_outer=cvpartition(y,'k',kFold);%

statistics = cell(CVP_outer.NumTestSets,16);% nroDatasets x nroModels to save the AUC values

%outer loop
for cv=1:CVP_outer.NumTestSets
 
    %training and test sets
    index=CVP_outer.training(cv);
    training=[{dataset{1}(index)} {dataset{2}(index,:)} {dataset{3}(index)} {dataset{4}(index)} {dataset{5}(index)} ];
    yT=training{4};
    kappa=training{5};
    
    index=CVP_outer.test(cv);
    test=[{dataset{1}(index)} {dataset{2}(index,:)} {dataset{3}(index)} {dataset{4}(index)} {dataset{5}(index)} ];
    
    yTest=test{4};

    %split the training set in training and validation sets
    CVP_inner=cvpartition(yT,'k',kFold);%
    
    [Q1,Q3,Q5,Q7,Q9]=findLambda(training{1});
    NN=length(training{1});
    grid_gamma=[1/Q1,1/Q3,1/Q5,1/Q7,1/Q9 ];
    grid_C=[1 1./(NN*[0.1:0.05:0.2])]; %no fraction of oultiers, 10% 15% and 20% of outliers
    
    for ModelOp=[1,6:16]
   % for ModelOp=[1:16] %for SMDDCCP (1), SMDDDA(2) SMDDDASN(3) OCSMM(4) SVDD(5);SMDDCCP (6-10) with k=0.90, .92, .94, .96, .98
        [cv ModelOp]
        %model selection with grid search
        [bestC,bestGamma, ~]=gridSearch(CVP_inner, training,ModelOp, grid_C, grid_gamma,3); %inner_loop for model selection
       
       
        %train the classifier with best hyperparameters
        %-----------------------------------------------
       [ kernelMatrices,X,Z,kappa] = getModelSetup( ModelOp,training,test,kappa, bestGamma,kernelOp);
        %PREDICT
        %--------
        ypred=prediction_MONQP(kernelMatrices, X, Z,bestC,kappa,kernelOp,ModelOp);
        %[ypred,yV]
        [~, ~,~,~,bestAUC,~]=Error_count(yTest,ypred);
        statistics{cv,ModelOp}=bestAUC;
        
    end
    
end
%modify this
save (['output/workspace' type num2str(option) '.mat'])
end

