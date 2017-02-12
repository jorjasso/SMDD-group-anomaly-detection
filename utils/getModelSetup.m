function [ kernelMatrices,X,Z,kappa] = getModelSetup( ModelOp,training,test,kappa, gamma,kernelOp)
%getModelSetup This procedure for a given model set the kernels "kernelMatrices", training
%and validation or (test sets) "X,Z" and kappa values "kappa".
%This procedure is very especific of the experiments, and it must be
%modified to include more models

% kernels for SMDDCPP, SMDD and OCSMM
%-------------------------------------------
if (ModelOp==1||ModelOp==2||ModelOp==4||ModelOp==6||ModelOp==7||ModelOp==8||ModelOp==9||ModelOp==10)
    X=training{1}; Z=test{1};
    [Krr,Kre,Kee,K_traceTraining,K_traceTest]=kernelPM(X, Z,kernelOp,gamma);
    kernelMatrices={Krr,Kre,Kee,K_traceTraining,K_traceTest};
    if (ModelOp==1)   kappa=ones(length(kappa),1); end
    if (ModelOp==6)   kappa=0.9*ones(length(kappa),1); end
    if (ModelOp==7)   kappa=0.8*ones(length(kappa),1); end
    if (ModelOp==8)   kappa=0.7*ones(length(kappa),1); end
    if (ModelOp==9)   kappa=0.6*ones(length(kappa),1); end
    if (ModelOp==10)  kappa=0.5*ones(length(kappa),1); end
    
    
end
% kernels for  SMDDCPP, SMDD with spherical normalization (option 4)
%-------------------------------------------
if (ModelOp==3 ||ModelOp==11||ModelOp==12||ModelOp==13||ModelOp==14||ModelOp==15||ModelOp==16)
    X=training{1}; Z=test{1};
    [Krr,Kre,Kee,K_traceTraining,K_traceTest]=kernelPM(X, Z,4,gamma);
    kernelMatrices={Krr,Kre,Kee,K_traceTraining,K_traceTest};
    % [Krr_SN,Kre_SN,Kee_SN,~]=kernelPM(X, Z,4,gamma);
    % [Krr,   Kre,   Kee,   K_traceTraining,K_traceTest]
    % kernelMatrices={Krr_SN,Kre_SN,Kee_SN};
    if (ModelOp==11)   kappa=ones(length(kappa),1); end
    if (ModelOp==12)   kappa=0.9*ones(length(kappa),1); end
    if (ModelOp==13)   kappa=0.8*ones(length(kappa),1); end
    if (ModelOp==14)   kappa=0.7*ones(length(kappa),1); end
    if (ModelOp==15)   kappa=0.6*ones(length(kappa),1); end
    if (ModelOp==16)  kappa=0.5*ones(length(kappa),1); end
end
%kernels for SVDD
%-------------------------------------------
if ModelOp==5
    X=training{2}; Z=test{2};
    Krr_SVDD=kernel(X,X,kernelOp,gamma);
    Kre_SVDD=kernel(X,Z,kernelOp,gamma);
    if kernelOp==3% RBF kernel
        Kee_SVDD=ones(size(Z,1),1);
    else
        Kee_SVDD=diag(kernel(Z,Z,kernelOp,gamma))';
    end
    kernelMatrices={Krr_SVDD,Kre_SVDD,Kee_SVDD};
end


end

