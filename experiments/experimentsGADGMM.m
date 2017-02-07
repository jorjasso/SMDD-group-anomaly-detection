function experimentsGADGMM(option,percentAnomalies,N)
%This function reproduce the experiments for the point-based group anomaly
%detection.
%option is a value in {1,2,3} depending on the type of anomalies.
%If option=4 then the anomalous groups are formed by combinig one third of
%each type {1,2,3} into one group. If option =5 the experiments the
%anomalies are distibution-based
%
%N is the size of the dataset.  (number of groups)
%
%percentAnomalies is the percent of anomalies %N

%%%%%%%%%%%%%%%%%%%%%%%%
%EXPERIMENT DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%
type='PointBasedGroupAnomalies';
kFold=10;
kernelOp=3
OpPrint=0
%%%%%%%%%%%%%%%%%%%%%%%%
%DATA
%%%%%%%%%%%%%%%%%%%%%%%%
switch option
    case 1 % point-based group anomalies first type
        [nonAnomalous, anomalous,~,~]=getPointBasedGroupAnomalies(N,OpPrint);
    case 2 % point-based group anomalies second type
        [nonAnomalous, ~,anomalous,~]=getPointBasedGroupAnomalies(N,OpPrint);
    case 3 % point-based group anomalies third type
        [nonAnomalous, ~,~,anomalous]=getPointBasedGroupAnomalies(N,OpPrint);
    case 4 % point-based group anomalies combined
        [nonAnomalous, A1,A2,A3]=getPointBasedGroupAnomalies(2*N,OpPrint);
        %the nonanomaous groups are formed by combining 1/3 o groups of
        %each A1, A2 and A3
        index=randperm(2*N,floor(2*N/3));
        T1={A1{1}(index),A1{2}(index,:),A1{3}(index)}
        T2={A2{1}(index),A2{2}(index,:),A2{3}(index)}
        T3={A3{1}(index),A3{2}(index,:),A3{3}(index)}
        T=[T1;T2;T3]
        anomalous{1}=[T{1,1} T{2,1} T{3,1}]
        anomalous{2}=[T{1,2}; T{2,2}; T{3,2}]
        anomalous{3}=[T{1,3} T{2,3} T{3,3}]
        
    case 5 % distribution based
        [nonAnomalous, anomalous]=getDistributionBasedData(N,OpPrint)
        
    otherwise
        disp('no valid option')
end

%percent anomalies
nAnomalies=floor(percentAnomalies/100*N);
nNonAnomalies=N-nAnomalies;

%clases
y=[ones(nNonAnomalies,1);-ones(nAnomalies,1)];

%kappa values for SMDD CCP
kappa=ones(length(y),1);


%index: selection nNonanomalies groups from nonAnomalous dataset and nAnomalies
%groups from anomalous dataset
index=randperm(N,nNonAnomalies)
NA={nonAnomalous{1}(index),nonAnomalous{2}(index,:),nonAnomalous{3}(index)}
index=randperm(N,nAnomalies)
A={anomalous{1}(index),anomalous{2}(index,:),anomalous{3}(index)}

%Combined dataset
tempData=[NA; A ]

dataset{1}=[tempData{1,1} tempData{2,1}]  %groups
dataset{2}=[tempData{1,2} ;tempData{2,2}] %means
dataset{3}=[tempData{1,3} tempData{2,3}]  %coovariance matrices
dataset{4}=y                              %classes
dataset{5}=kappa                          %kappa values

%%%%%%%%%%%%%%%%%%%%%%%%
%Nested cross validation
%%%%%%%%%%%%%%%%%%%%%%%%
CVP_outer=cvpartition(y,'k',kFold);%

statistics = cell(CVP_outer.NumTestSets,16);% nroDatasets x nroModels to save the AUC values

%outer loop
for cv=1:CVP_outer.NumTestSets
    cv
    %training and test sets
    index=CVP_outer.training(cv)
    training=[{dataset{1}(index)} {dataset{2}(index,:)} {dataset{3}(index)} {dataset{4}(index)} {dataset{5}(index)} ]
    yT=training{4};
    kappa=training{5};
    
    index=CVP_outer.test(cv)
    test=[{dataset{1}(index)} {dataset{2}(index,:)} {dataset{3}(index)} {dataset{4}(index)} {dataset{5}(index)} ]
    
    yTest=test{4};

    
    
    %split the training set in training and validation sets
    CVP_inner=cvpartition(yT,'k',kFold);%
    
    [Q1,Q3,Q5,Q7,Q9]=findLambda(training{1})
    grid_gamma=[1/Q1,1/Q3,1/Q5,1/Q7,1/Q9 ]
    grid_C=[1 1./(nTraining*[0.1:0.1:0.2])]
    
    %temp borrar
    grid_gamma=[1/Q1,1/Q3]
    grid_C=[1 1./(nTraining*[0.1])]
    %
    
    
    for ModelOp=[1:16] %for SMDDCCP (1), SMDDDA(2) SMDDDASN(3) OCSMM(4) SVDD(5);SMDDCCP (6-10) with k=0.90, .92, .94, .96, .98
        ModelOp
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

