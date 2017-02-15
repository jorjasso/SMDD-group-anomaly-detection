% jorge.jorjasso@gmail.com,
% modified version of gridSearch. m in TSK-kernels-Low-quality
%  at https://github.com/jorjasso/TSK-kernels-Low-quality/blob/master/gridSearch.m
% date 09 set 2015
function [bestC,bestGamma, bestAUC]=gridSearch(CVP, training,ModelOp, grid_C, grid_gamma,kernelOp)

kFold=CVP.NumTestSets;
history_AUC=zeros(1,kFold);
grid_AUC=zeros(21,21);

j=1;

%for testing
%grid_gamma=1
%grid_C=1
%kFold=1
%

for gamma=grid_gamma
    gamma
    i=1;
    
    for C=grid_C
        
        for k = 1:kFold
            %[gamma, C, k]
            indTraining=CVP.training(k);
            indValidation=CVP.test(k);
            
            %subset of the data to be used for training and validation
            %---------------------------------------------------------
            trainCV=[{training{1}(indTraining)} {training{2}(indTraining,:)} {training{3}(indTraining)} {training{4}(indTraining)} ];
            validationCV=[{training{1}(indValidation)} {training{2}(indValidation,:)} {training{3}(indValidation)} {training{4}(indValidation)} ];
            
            yT=trainCV{4};
            yV=validationCV{4};
            
            %As one-class classifiers does not have labels we train the
            %classifiers using only non-anomalous points, and validate with
            %anomalous and anomalies. By optimizing the hyper-paramters for
            %the AUC we think that we will estimate the right description
            %for the data
            
            %Redefining the training set: Using only the non-anomalous points for training the
            %one-class classifieres
            index=(yT==1);
            trainCV=[{training{1}(index)} {training{2}(index,:)} {training{3}(index)} {training{4}(index)} ];
            
            %yT contain only points with label = 1
            
           [ kernelMatrices,X,Z,kappa] = getModelSetup( ModelOp,trainCV,validationCV, gamma,kernelOp);

            %predict
            
            ypred=prediction_MONQP(kernelMatrices, X, Z,C,kappa,kernelOp,ModelOp);
           
            [~, ~,~,~,AUC,~]=Error_count(yV,ypred);
            
            history_AUC(k)=AUC;
        end
        %grid containing the mean of cross validated AUC values
        grid_AUC(i,j)=mean(history_AUC);
        
        i=i+1;
    end
    j=j+1;
end

[index_C,index_gamma]=bestParammeters(grid_AUC);

bestC     = grid_C(index_C);
bestGamma = grid_gamma(index_gamma);
bestAUC=grid_AUC(index_C,index_gamma);
end



%----------------------Search the best parammeteres---return the indexes---------------------------------------
function [f,c]=bestParammeters(grid)
[fil,col]=find(grid==max(grid(:)));
f=fil(1);
c=col(1);
end