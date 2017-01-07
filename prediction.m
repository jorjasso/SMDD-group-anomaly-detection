function statistics=prediction(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
% prediction on the test set for several models
S=training{1};muTrain=training{2};Sigma=training{3};
STest=test{1}; muTest=test{2};SigmaTest=test{3};

if ModelOp==1 [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam); end
if ModelOp==2 [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,kernelOp,kernelParam) ; end
if ModelOp==3 [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,4,kernelParam) ; end  %spherical normalization option 4
if ModelOp==4 [rho,alph1,Io,I1,I2,I3]=OCSMM(S,muTrain,Sigma,C,kernelOp,kernelParam) ; end
if ModelOp==5 [R,c,cNorm,alph1,Io,I1,I2,I3]=SVDD(muTrain,C,kernelOp,kernelParam);  end

if (kernelOp~= 4 & ModelOp~=5)  %if is no gaussian kernel with spherical normalization or SVDD
    [~,Kgrid]=kernelPM(S(Io),muTrain(Io,:),Sigma(Io),STest,muTest,SigmaTest,kernelOp,kernelParam);
    [K_trace,KTest]=kernelPM(STest,muTest,SigmaTest,STest,muTest,SigmaTest,kernelOp,kernelParam);
end

% 1 anomalous, -1 normal
e=ones(length(STest),1);
if ModelOp==1 ypred= diag(KTest)-(2*Kgrid'*alph1(Io))/sum(alph1(Io))+cNorm.*e+K_trace-R^2*e ; end
if ModelOp==2 ypred= diag(KTest)-(2*Kgrid'*alph1(Io))+cNorm.*e-R^2*e; end
if ModelOp==3
    [~,Kgrid]=kernelPM(S(Io),muTrain(Io,:),Sigma(Io),STest,muTest,SigmaTest,4,kernelParam); %spherical normalization option 4
    [K_trace,KTest]=kernelPM(STest,muTest,SigmaTest,STest,muTest,SigmaTest,4,kernelParam); %spherical normalization option 4
    ypred= diag(KTest)-(2*Kgrid'*alph1(Io))+cNorm.*e-R^2*e; end
if ModelOp==4 ypred=rho-Kgrid'*alph1(Io); end
if ModelOp==5
    Kgrid=kernel(muTrain(Io,:),muTest,kernelOp,kernelParam);
    KTest=kernel(muTest,muTest,kernelOp,kernelParam);
    ypred=diag(KTest)-(2*Kgrid'*alph1(Io))+cNorm.*e-R^2*e; end  
[Err_Rate, Err_RateA,Err_RateN,~,AUC,stat]=Error_count(yt,ypred);
statistics=[Err_Rate AUC  Err_RateA Err_RateN length(Io) length(I1) length(I2) C stat(12) stat(13)];
