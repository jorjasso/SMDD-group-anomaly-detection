run   /home/jorjasso/cvx/cvx_startup.m

nTraining=50;
nValidation=30;
nTest=30;
kernelOp=3;
kappa=ones(nTraining,1);
lambdaValues=[];
yv=[ones(20,1);-ones(10,1)]; %there are 20 anomalous groups and 10 normal groups for validation
yt=[ones(20,1);-ones(10,1)]; %there are 20 anomalous groups and 10 normal groups for test
for i=1:200
    
    [training, validation, test]=generateAnomalusGroupDetectionDataMultimodalII(nTraining,nValidation,nTest,0)
    %kernel parameter
    S=training{1};
    kernelParam=findLambda(S);
    lambdaValues=[lambdaValues kernelParam];
    
    %SMDDCPP model
    %----------------------------------------------------------------------   
    
    for ModelOp=[1:5] %for SMDDCCP, SMDDDA SMDDDASN OCSMM SVDD
        %Model selection based on AUC over C
        %----------------------
        MetricOp=1;
        C_AUC=findC(training,validation,kappa,kernelOp,kernelParam,ModelOp,MetricOp,yv);
        
        %TEST, goal metric AUC
        %--------------
        statAUC=prediction(training, test, kernelParam,C_AUC,kappa,kernelOp,ModelOp,yt);
        statistics_AUC{i,ModelOp}=statAUC;
        
        %Model selection based on ACC over C
        %------------------------------------
        MetricOp=2;
        C_ACC=findC(training,validation,kappa,kernelOp,kernelParam,ModelOp,MetricOp,yv);
        
        %TEST, goal metric ACC
        %--------------
        statACC=prediction(training, test, kernelParam,C_ACC,kappa,kernelOp,ModelOp,yt);
        statistics_ACC{i,ModelOp}=statACC;
    end        
end

save experimentosGroupAnomaluGMMMultimodal.mat
