function C=findC(training,validation,kappa,kernelOp,kernelParam,ModelOp,MetricOp,yt)
% function C=findC(training,kappa,kernelOp,kernelParam,ModelOp), performrs model
% selecion over C, training is a cell of the form training[replicates, mean(replicates) cov(replicates)]
% ModelOp is the model option: ModelOp=1 SMDDCCP, ModelOp=2 SMDDDA, ModelOp=3 SMDDDA Spherical
% normalization, ModelOp=4 OCSMM, ModelOp=5 SVDD
% MetricOp =1 model selection based on AUC,
% MetricOp =2 model selection based on Accuracy,

statistics=[];
for C=2.^[-5:1:0]
    stat=prediction_MONQP(training, validation, kernelParam,C,kappa,kernelOp,ModelOp,yt);
    statistics=[statistics;stat];    
end

if MetricOp==1
    ind=find(statistics(:,2)==max(statistics(:,2))); % for several max indices
    ind2=find(statistics(ind,1)==min(statistics(ind,1)));% chose the AUC with minimum error rate, (ind2)
    C=max(statistics(ind(ind2),8));% for several ind2 indices, chose  the maximum C
    
end

if MetricOp==2
    %chose best model based on ACC
    %----------------
    ind=find(statistics(:,1)==min(statistics(:,1))); % for several min indices (search the minimum error)
    ind2=find(statistics(ind,2)==max(statistics(ind,2)));% chose the AUCC with maximum AUC, (ind2)
    C=max(statistics(ind(ind2),8));% for several ind2 indices, chose  the maximum C
end

