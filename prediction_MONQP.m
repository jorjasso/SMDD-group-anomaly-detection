function statistics=prediction_MONQP(kernelMatrices,S, STest,C,kappa,kernelOp,ModelOp,yt)
% prediction  for several models
Krr=kernelMatrices{1};
Kre=kernelMatrices{2};
Kee=kernelMatrices{3};
if ModelOp==1
    K_traceTraining=kernelMatrices{4};
    K_traceTest=kernelMatrices{5};
end

%training
%---------
if ModelOp==1 [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(Krr,K_traceTraining',S,kappa,C,kernelOp); end
if (ModelOp==2|| ModelOp==3) [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA_MONQP(Krr,S,C,kernelOp) ; end
%if ModelOp==3|| [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA_MONQP(Krr,S,C,kernelOp,kernelParam) ; end  %spherical normalization option 4
if ModelOp==4 [rho,alph1,Io,I1,I2,I3]=OCSMM_MONQP(Krr,S,C,kernelOp) ; end
if ModelOp==5 [R,c,cNorm,alph1,Io,I1,I2,I3]=SVDD_MONQP(Krr,S,C,kernelOp);  end

%test ( 1 anomalous, -1 normal)
%-------
e=ones(length(STest),1);

if ModelOp==1 ypred= Kee'-(2*Kre(Io,:)'*alph1(Io))/sum(alph1(Io))+cNorm.*e+K_traceTest'-R^2*e ; end
if (ModelOp==2|| ModelOp==3) ypred= Kee'-(2*Kre(Io,:)'*alph1(Io))+cNorm.*e-R^2*e; end
%if ModelOp==3 ypred= Kee-(2*Kre'*alph1(Io))+cNorm.*e-R^2*e; end
if ModelOp==4 ypred=rho-Kre(Io,:)'*alph1(Io); end
if ModelOp==5 ypred=Kee-(2*Kre(Io,:)'*alph1(Io))+cNorm.*e-R^2*e; end  

[Err_Rate, Err_RateA,Err_RateN,Borrar,AUC,stat]=Error_count(yt,ypred);
statistics=[Err_Rate AUC  Err_RateA Err_RateN length(Io) length(I1) length(I2) C stat(12) stat(13)];

