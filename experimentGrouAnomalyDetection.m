function [Err_Rate ConfMat Io I1 I2] = experimentGrouAnomalyDetection(kernelParam, C,kappa)

muGM = [-1.7 1; 1.7 -1; 0 2];
I=0.2*eye(2);
%sigma = cat(3,0.2*I,0.2*I,0.2*I);
sigmaGM = I;
p=ones(1,3)/3;%weigths for normal groups
obj = gmdistribution(muGM,sigmaGM,p);

%generate 50 groups with number of points per group N~poisson(100)

for i=1:50
    n=poissrnd(100);
    S{i} = random(obj,n);        
    mu(i,:)=mean(S{i}); %it is not used actually
    Sigma{i}=cov(S{i}); %it is not used actually        
end
%actually mu and Sigma are not used, because kernelPM.m and computeKernel
%use only information of replicates S{i} to compute the kernel on
%probability measures and the trace, however mu and Sigma could be used in
%anohter kernel settings on pm, see computekernel.m case 5 and so on




%group of point anomalies
i=1;
n=poissrnd(100);
STest{i} = mvnrnd([0,0],I,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});

%group anomaly for points individually normal
i=i+1;
p=[0.85,0.08,0.07];%weigths for anomalous groups
obj = gmdistribution(muGM,sigmaGM,p);
n=poissrnd(100);
STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});

i=i+1;
p=[0.04,0.48,0.48];%weigths for anomalous groups
obj = gmdistribution(muGM,sigmaGM,p);
n=poissrnd(100);
STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});


n=50;
kernelOp=3; %gaussian kernel



% make model selection good parameters
%kernelParam=2^(-5); C=2^(1);
%kernelParam=2^(-4); C=2^(1);
%kernelParam=2^(-3); C=2^(1);
%kernelParam=2^(4); C=2^(-1);
%{C=2^{-5} ,kernelParam=2^(0)} {C=2^(-4),}kernelParam=2^(3)};
% {C=2^(-3);kernelParam=2^(4);}  {C=2^(-2);kernelParam=2^(4);}
% {C=2^(-1);kernelParam=2^(4);}
[R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,mu,Sigma,kappa,C,kernelOp,kernelParam);
[~,Kgrid]=kernelPM(S(Io),mu(Io),Sigma(Io),STest,muTest,SigmaTest,kernelOp,kernelParam);
[K_trace,KTest]=kernelPM(STest,muTest,SigmaTest,STest,muTest,SigmaTest,kernelOp,kernelParam);
% -1 anomalous, 1 normal
e=ones(3,1);
score_p=R^2*e -diag(KTest)+(2*Kgrid'*alph1(Io))/sum(alph1(Io))-cNorm.*e-K_trace; %poner solo con los sv
predict_label=sign(score_p);
[Err_Rate, ConfMat]=Error_count([-1;-1;-1],predict_label)



[X,Y,T,AUC] = perfcurve([-1;-1;-1], score_p, -1);
plot(X,Y)

%ROC analysis
testpred=score_p;
testt=[-1;-1;-1];
th_vals = [min(testpred):0.01:max(testpred)+0.01];
sens = [];
spec = [];
for i = 1:length(th_vals)
    b_pred = testpred>=th_vals(i);
    TP = sum(b_pred==1 & testt == 1);
    FP = sum(b_pred==1 & testt == -1);
    TN = sum(b_pred==0 & testt == -1);
    FN = sum(b_pred==0 & testt == 1);
    sens(i) = TP/(TP+FN);
    spec(i) = TN/(TN+FP);
end

figure(1);hold off
cspec = 1-spec;
cspec = cspec(end:-1:1);
sens = sens(end:-1:1);
plot(cspec,sens,'k')
AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));



