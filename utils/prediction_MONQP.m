function ypred=prediction_MONQP(kernelMatrices,S, STest,C,kappa,kernelOp,ModelOp)

% prediction  for several models
Krr=kernelMatrices{1};
Kre=kernelMatrices{2};
Kee=kernelMatrices{3};
if (ModelOp==1 ||ModelOp==6 ||ModelOp==7||ModelOp==8||ModelOp==9||ModelOp==10||ModelOp==11||ModelOp==12||ModelOp==13||ModelOp==14||ModelOp==15||ModelOp==16)

    K_traceTraining=kernelMatrices{4};
    K_traceTest=kernelMatrices{5};
end

%training
%---------
if (ModelOp==1 ||ModelOp==6 ||ModelOp==7||ModelOp==8||ModelOp==9||ModelOp==10||ModelOp==11||ModelOp==12||ModelOp==13||ModelOp==14||ModelOp==15||ModelOp==16)
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(Krr,K_traceTraining',S,kappa,C,kernelOp); 
end
if (ModelOp==2|| ModelOp==3) 
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA_MONQP(Krr,S,C,kernelOp) ; 
end
if ModelOp==4 
    [rho,alph1,Io,I1,I2,I3]=OCSMM_MONQP(Krr,S,C,kernelOp) ; 
end
if ModelOp==5 
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SVDD_MONQP(Krr,S,C,kernelOp);
end

%test ( 1 anomalous, -1 normal)
%-------
e=ones(length(STest),1);

if (ModelOp==1 ||ModelOp==6 ||ModelOp==7||ModelOp==8||ModelOp==9||ModelOp==10||ModelOp==11||ModelOp==12||ModelOp==13||ModelOp==14||ModelOp==15||ModelOp==16)
     ypred= Kee'-(2*Kre(Io,:)'*alph1(Io))/sum(alph1(Io))+cNorm.*e+K_traceTest'-R^2*e ; end
if (ModelOp==2|| ModelOp==3) ypred= Kee'-(2*Kre(Io,:)'*alph1(Io))+cNorm.*e-R^2*e; end
if ModelOp==4 ypred=rho-Kre(Io,:)'*alph1(Io); end
if ModelOp==5 ypred=Kee-(2*Kre(Io,:)'*alph1(Io))+cNorm.*e-R^2*e; end  


