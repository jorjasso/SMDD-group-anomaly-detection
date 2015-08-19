%% kernel version SVDDPM
%  modified 12 mar 2014
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
clc;clear;
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)))
% artificial data
n=10    ;d=2;L=5;
[S,muTrain,Sigma]=getData(n,d,L);

C=1;
kappa=0.8 + (1-0.8).*rand(n,1);
kappa=ones(n,1)*1;

%kernelOp=1; kernelParam=[]%lineal
%kernelOp=2; kernelParam=[0,3]; %polinomial
kernelOp=3; kernelParam=[2^(1),0]; %gausiaan


[R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)

%-------------------------------------------------------------

%Experiment plot circle for USVDD with linear kernel
%  X=muTrain;
%
%  plotCircle(X,Sigma,R,c,0,'-r',10)
%  trSigma=ones(n,1)
%  for i=1:n
%      trSigma(i)=trace(Sigma{i});
%  end
%
% % %EXPERIMENT
% % %USVDD and USVDDKernel gives the same solutions (R, R1) and (c
% % %c1) for the case of linear kernel: test ok
%  [R1,c1,alph11,Io1,Isurface,Ierrors,IinSphere]=USVDD(X,trSigma,kappa,C)
%
%  plotCircle(X,Sigma,R1,c1,0,'-r',10)

%% Figure, different values for gaussian kernel param, same C
kernelOp=3; kernelParam=[2^(1),0]; %gausiaan
C=1;
for par=2.^[-10:1:10]
    kernelParam=[par,0]
    [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)
end

%% Figure different values for gaussian kernel param, vs kappa, same C
% 1)
C=1;
kappa=ones(n,1)*1;
j=1;
for kappa=ones(n,1)*[0.7 0.8 0.9 1]
    for par=2.^[-3:1:3]
        kernelParam=[par,0]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)
        
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
end

%  the same (1) but plot only for some values of kappa and gammaRBF
%2)
j=1;
for kappa=ones(n,1)*[0.8 1]
    for par=2.^[-1 0 2]
        kernelParam=[par,0]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)
        
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
end

% all kappa values constant menos dos ambos con kappa=0.8 y kappa=0.9

for kappa=ones(n,1)*[1]
    for par=2.^[ 0 ]
        for val=[1 0.9 0.8]
            kappa(7)=val;
            kappa(8)=val;
            kernelParam=[par,0]
            [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
            plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)
            
            saveas(gcf,['copiar/fig' int2str(j)],'fig');
            print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
            close(gcf);
            j=j+1;
        end
    end
end

clear;
clc;
load data.mat

%the same as 2) but with a large covariance matrix,
% result: the trace in the RKHS is a small value, is for that is no more
% difference for kappa=1 and SVDD
for i=1:n
    Sigma{i}=Sigma{i}*100;
end


j=1;
for kappa=ones(n,1)*[1]
    for par=2.^[-1 0 2]
        kernelParam=[par,0]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
        %            plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam,0)
        %            figure
        [R1,c1,cNorm1,alph,Io,I1,I2,I3]=SMDDPMDA(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam)
        %            plotContourKernelPM(S,muTrain,Sigma,R1,cNorm1,alph,Io,I1,kernelOp,kernelParam,1)
        [R R1 abs(R-R1)]
        %saveas(gcf,['copiar/fig' int2str(j)],'fig');
        %print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        %close(gcf);
        %j=j+1;
    end
end


%% Experiments for SMDDDA vs SMDDCCP
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
% 1) linear kernel
clear;
clc;
kernelOp=1
n=10; d=2;L=5;
kappa=ones(n,1)
C=1;
[S,muTrain,Sigma]=getData(n,d,L,0.5);



clear;
load datosSMDDDAvsSMDDDCPP
j=1;
matrR=[];

for par=2.^[ 0] % it has no efecto because is the linear kernel
    kernelParam=[par,0]
    [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam,0)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    
    [R1,c1,cNorm1,alph,Io,I1,I2,I3]=SMDDPMDA(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam)
    plotContourKernelPM(S,muTrain,Sigma,R1,cNorm1,alph,Io,I1,kernelOp,kernelParam,1)
    matrR=[matrR; R R1 abs(R-R1)];
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
end


% 2) gaussian kernel


clear;
clc;
kernelOp=3
n=10; d=2;L=5;
kappa=ones(n,1)
C=1;
[S,muTrain,Sigma]=getData(n,d,L,0.5);

j=1;
matrR=[];

for par=2.^[ -1 0 2]
    for kappa=ones(n,1)*[0.8 0.9 1]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0)
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    
    
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,kappa,C,kernelOp,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,1)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
end

j=1;
for par=2.^[ -4:0.3:-1]
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,4,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,1)
    %     saveas(gcf,['copiar/fig' int2str(j)],'fig');
    %     print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    %     close(gcf);
    %     j=j+1;
end
%experimentar con  par=2^(-3), parece que con el kernel lineal ok, pero
%con el kernel gaussiano el radio del modelo que no considera el trace
%crece demasiado
% experimentar con los otros valores
j=1;
for par=2.^[ -15:1:2]
    for kappa=ones(n,1)*[ 0.7 0.8 0.9 1]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
        
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0)
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    j=j+1;
end

%-------------------------------------------------------------------------------------------
%experimentos con la misma matrix de convarianza para todos los datos,
%la idea es testar no usando dirac distribution sino usando
%distribuciones con la misma covarianza y ver que pasa con los dos
for i=1:n
    Sigma{i}=Sigma{1};
end

j=1;
for par=2.^[ -7:1:7]
    for kappa=ones(n,1)*[0.8 0.9 1]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0)
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    %kernelOp=4, spherical normalization
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,4,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,4,par,1)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    [rho,alph1,Io]=OCSMM(S,muTrain,Sigma,C,kernelOp,par)
    plotContourKernelPM(S,muTrain,Sigma,rho,-1,alph1,Io,I1,kernelOp,par,2)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
end

run  run /home/jorjasso/cvx/cvx_startup.m
clc
clear
load    SMDDDAvsSMDDCCPGaussian
save SMDDDAvsSMDDCCPGaussian.mat


%% experiments in copiar 1run  run /home/jorjasso/cvx/cvx_startup.m
run   /home/jorjasso/cvx/cvx_startup.m
clc
clear
load    SMDDDAvsSMDDCCPGaussian
j=1;
for par=2.^[ -2:1:7]
    for kappa=ones(n,1)*[0.8 0.9 1]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0)
        saveas(gcf,['copiar1/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar1/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    %kernelOp=4, spherical normalization
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,4,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,4,par,1)
    saveas(gcf,['copiar1/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar1/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    [rho,alph1,Io]=OCSMM(S,muTrain,Sigma,C,kernelOp,par)
    plotContourKernelPM(S,muTrain,Sigma,rho,-1,alph1,Io,I1,kernelOp,par,2)
    saveas(gcf,['copiar1/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar1/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
end


%%  figures different values for kappa vs C, same parm of kern RBF
n=20;d=2;L=5;
[S,muTrain,Sigma]=getData(n,d,L,0.3);
kappa=ones(n,1)*1;
j=1;
kernelParam=[2^0 0];
nSV=ones(4,6)*-1
nErr0=ones(4,6)*-1
r=1;s=1;
for kappa=ones(n,1)*[0.7 0.8 0.9 1]
    for C = 2.^[-3:0.2:-2]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=USVDDKernel(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,kernelParam)
        
        nSV(r,s)=length(I1);
        nErr(r,s)=length(I2);
        s=s+1;
        
        
        ['\it{kappa' int2str(kappa(1)) ' vs C} C =' int2str(C)]
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    r=r+1;
    s=1;
end

%% modify the data set of replicates to make a nice plot april 3 2014
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSMDDDAvsSMDDDCPP
%modify the data set of replicates to make a nice plot
%------------------------------------------
L=15;
Sigma{2}=Sigma{2}*5.7;
muTrain(2,1)=-0.5
S{2} = mvnrnd(muTrain(2,:),Sigma{2},L);
muTrain(2,:)=mean(S{2});             % estimated mean
Sigma{2}=cov(S{2});

muTrain(3,2)=3.5;
i=3;
S{i} = mvnrnd(muTrain(i,:),Sigma{i},L);
muTrain(i,:)=mean(S{i});             % estimated mean
Sigma{i}=cov(S{i});


muTrain(4,:)=[0.7 1.5];
i=4;
S{i} = mvnrnd(muTrain(i,:),Sigma{i},L);
muTrain(i,:)=mean(S{i});             % estimated mean
Sigma{i}=cov(S{i});


muTrain(6,2)=5;
i=6;
S{i} = mvnrnd(muTrain(i,:),Sigma{i},L);
muTrain(i,:)=mean(S{i});             % estimated mean
Sigma{i}=cov(S{i});


muTrain(7,1)=0.5;
i=7;
S{i} = mvnrnd(muTrain(i,:),Sigma{i},L);
muTrain(i,:)=mean(S{i});             % estimated mean
Sigma{i}=cov(S{i});


muTrain(9,:)=[1.9 2.7 ];
i=9;
S{i} = mvnrnd(muTrain(i,:),Sigma{i},L);
muTrain(i,:)=mean(S{i});             % estimated mean
Sigma{i}=cov(S{i});

save datosSMDDDAvsSMDDDCPPModificados.mat
%------------------------------------------
%% SMDDA vs SMDDCCP with linear kernel april 3 2014
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSMDDDAvsSMDDDCPPModificados

kernelOp=1
kappa=ones(n,1)
C=1;

%covariance for dirac test measures
[~,d]=size(S{1});
SigmaT=0.001 + (0.01-0.001).*rand(d,d);  %initial Sigmas
SigmaT=(SigmaT+SigmaT'+eye(d));
SigmaTestGrid=SigmaT*0.0000009; %small covariance matrix to aproximate Dirac measures

%SMDDCCP
[R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0,SigmaTestGrid)
saveas(gcf,['copiar/figSMDDCCP'],'fig');
print(gcf, '-dpng', '-r0', ['copiar/figSMDDCCP.png']);
close(gcf);

%SMDDDA
[R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,kernelOp,par)
plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,1,SigmaTestGrid)
saveas(gcf,['copiar/figSMDDDA'],'fig');
print(gcf, '-dpng', '-r0', ['copiar/figSMDDDA.png']);
close(gcf);


%OCSVM
[rho,alph1,Io,I1]=OCSMM(S,muTrain,Sigma,C,kernelOp,par)
plotContourKernelPM(S,muTrain,Sigma,rho,-1,alph1,Io,Io,kernelOp,par,2,SigmaTestGrid)
saveas(gcf,['copiar/figOCSVM'],'fig');
print(gcf, '-dpng', '-r0', ['copiar/figOCSVM.png']);
close(gcf);

%% Experiment: Role of kappa values  10 april
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSMDDDAvsSMDDDCPPModificados
load datosSameSigma
kernelOp=3, par=2^(0);
kappa=ones(n,1)
C=1;
%training set Same Sigma, Test set grid of pm with the same sigma of the
%training set
SigmaTestGrid=Sigma{1};
j=1
for i=[1 0.9 0.8 0.7]
    kappa=i*ones(n,1)
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0,SigmaTestGrid)
    saveas(gcf,['copiar/figSameSigmaTestGrid' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/figSameSigmaTestGrid' int2str(j) '.png']);
    close(gcf);
    j=j+1;
end

%training set Same Sigma, Test set grid of Dirac pm
[~,d]=size(S{1});
SigmaT=0.001 + (0.01-0.001).*rand(d,d);  %initial Sigmas
SigmaT=(SigmaT+SigmaT'+eye(d));
SigmaTestGrid=SigmaT*0.0000009; %small covariance matrix to aproximate Dirac measures

j=1
for i=[1 0.9 0.8 0.7]
    kappa=i*ones(n,1)
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0,SigmaTestGrid)
    saveas(gcf,['copiar/figDiracMeasureTest' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/figDiracMeasureTest' int2str(j) '.png']);
    close(gcf);
    j=j+1;
end

%conclusion: use Dirac distribution to vizualice.
load datosSameSigmaIV

kappa=ones(n,1)
j=1
kernelOp=3, par=2^(-1);
for val=[0.8 0.9 1]
    kappa=val*ones(n,1);
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0,SigmaTestGrid)
    saveas(gcf,['copiar/figDiracMeasureTestkappaSigma125G' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/figDiracMeasureTestkappaSigma125G' int2str(j) '.png']);
    close(gcf);
    j=j+1;
end

kernelOp=3, par=2^(-1);
kappa=ones(n,1)
j=1
for val=[0.8 0.9 1]
    kappa(2)=val;
    kappa(1)=val;
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0,SigmaTestGrid)
    saveas(gcf,['copiar/figDiracMeasureTestkappaSigma125Diff' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/figDiracMeasureTestkappaSigma125Diff' int2str(j) '.png']);
    close(gcf);
    j=j+1;
end




%% SMDDA vs SMDDCCP with gaussian kernel april 3 2014
%modificado 07 april: SMDDCCP vs SMDDQP vs SMDDDA vs SMDDDA without
%normalization using a gaussian mixture model.
run /home/jorjasso/Downloads/cvx/cvx_startup.m

%---------------------------------------------------------------------
% First experiment: a gaussian mixture model with the same p_i, same
% covariance matrices
load datosSameSigma

kernelOp=3;
C=1;
j=1;
SigmaTestGrid=0.000001*[1 0;0 1];
for par=2.^[ -10:1:7]
    %SMDDCCP
    ypred=0;
    for kappa=ones(n,1)*[0.8 0.9 1]
        [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
        %[R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,1/sigma3);
        plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,ypred,SigmaTestGrid)
        saveas(gcf,['copiar/fig' int2str(j)],'fig');
        print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
        close(gcf);
        j=j+1;
    end
    
    %SMDDQP
    ypred=0;
    [rho,c,cNorm,alph,Io,I1,I2,I3]=SMDDCCPASOCSVM(S,muTrain,Sigma,kappa,C,kernelOp,par)
    [K_trace,~]=kernelPM(S(Io),muTrain(Io),Sigma(Io),S(Io),muTrain(Io),Sigma(Io),kernelOp,par);
    R=sqrt((K_trace(1)+cNorm-2*rho)/kappa(1));
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,ypred,SigmaTestGrid)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    
    %SMDDA
    ypred=1;
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,kernelOp,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,ypred,SigmaTestGrid)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    %kernelOp=4, spherical normalization
    %SMDDA with spherical normalization
    ypred=1;
    [R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDDA(S,muTrain,Sigma,C,4,par)
    plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,4,par,ypred,SigmaTestGrid)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
    %OCSVM
    ypred=2;
    [rho,alph1,Io,I1]=OCSMM(S,muTrain,Sigma,C,kernelOp,par)
    plotContourKernelPM(S,muTrain,Sigma,rho,-1,alph1,Io,I1,kernelOp,par,ypred,SigmaTestGrid)
    saveas(gcf,['copiar/fig' int2str(j)],'fig');
    print(gcf, '-dpng', '-r0', ['copiar/fig' int2str(j) '.png']);
    close(gcf);
    j=j+1;
    
end
%-----------------------------------------------------------------------

load datosSMDDDAvsSMDDDCPPModificados


%% Experiment over a GMM dataset
run   /home/jorjasso/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
addpath ./SVM-KM/


%Unimodal GMM
%------------
experimentGroupAnomalyDetectionGMM(0,200)

%Multimodal GMM
%--------------
experimentGroupAnomalyDetectionGMM(1,200)


% sloan resultados,
% los resultados estan dentro de la carpeta SVDDUNcertainData del dropbox
% me gustan:
% C=4;lambda=2; (load experimentosSDSS0Sigma) (prbar con lamnda 3 y lambda 4)
%matrixC = cell2mat(statistics(:,:,C,lambda)); %[epoch ModelOp C]
%------------------------------

% load experimentosSDSS1Sigma
% C=3;lambda=3;

% load experimentosSDSS5Sigma
% C=1;lambda=3; 1 4, 32, 33

% load experimentosSDSS10Sigma (1, 2) (1 3)

% experimentosSDSS100Sigma (1,4)
%---------------------------------
% GMM en drpbox...SVDUNCERTAINATA/ analizar los .mat :
% experimentosGroupAnomaluGMMUnimodal: tiene los experimentos del
% distribution based con escala [1 0; 1 0]
% experimentosGroupAnomaluGMMAnother: tiene los experimentos distribution
% based segun Xion
%---------------------------------

% escience0 /anomalousScenes  experiments ofr anomalous scences

% para lambda=9 (no lambda=1/sigma) funciona mejor para el multimodal, para
% el unimodal lambda y 1/lambda no mejoran, experimentar con los quantiles

% experimentos en escience0 /SMDD/anomalous/ images 10 iterations
kernelOp=5
statistics{epoch,ModelOp,j,k}=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt);

%example of files
%function d = L2_distance(a,b,df)
% L2_DISTANCE - computes Euclidean distance matrix
%
% E = L2_distance(A,B)
%
%    A - (DxM) matrix
%    B - (DxN) matrix
%    df = 1, force diagonals to be zero; 0 (default), do not force
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description :
%    This fully vectorized (VERY FAST!) m-file computes the
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example :
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Wed Oct 20 08:58:08 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

% Copyright notice: You are free to modify, extend and distribute
%    this code granted that the author of the original code is
%    mentioned as the original author of the code.

% Fixed by JBT (3/18/00) to work for 1-dimensional vectors
% and to warn for imaginary numbers.  Also ensures that
% output is all real, and allows the option of forcing diagonals to
% be zero.  

% trabalho con imagenes: escribir sobre low rank y codificar para anolaous
% scene detection
% Works as advertised! To train a linear SVM in the Nystrom feature-space, replace the computation of G and Ktilde at the end of INys() with:
% 
% M = Ve(:,pidx) * inVa;
% Mdata = E * M;
% 
% then train SVM with Mdata (which has all training samples in the new feature space). To compute the same features for a test vector z, use:
% 
% Mz = exp(-sqdist(z', center')/kernel.para) * M;
% 
% and classify with the same linear SVM.



%% SLOAN complete 500 dimensional data set
% experimentos for SDSS dataset, 
% are considered 200 runs, (epochs) to get statistics, 
% anomalous groups are generated from two different ways:
%     experiementSloan(0)= select randomly 50 clusters from the 505 normal
%     clusters,  to form the test set, additionally, group anomalies are
%     injected by generating 50 groups, each group with random galaxies.
%     then the test set has 100 groups.
% 
%     experiementSloan(factor), where factor is a integer number different of
%     zero, select 50 normal group from the total fo clusters, and group
%     anmalies are injected from a 4-dimensional GMM with three components, with
%     parameters: p=1/3; the means are randomly selected using information of the dataset
%     (the means of each group are selected using poisson(10) random points
%     of the dataset see experiementSloan.m lines 87-95, ) and the covariance of GMM is factor*sigmaMean, 
%     sigmaMean, is the mean of the covariance matrix of the clusters.
%      
%    information of statistics:
%       statistics=[Err_Rataddpath ./anomalousScenes/SceneClass13/MITcoast
addpath ./anomalousScenes/SceneClass13/MITforest
addpath ./Improved_Nystrom_Method

run   /home/jorjasso/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
addpath ./SVM-KM/
e AUC  Err_RateA Err_RateN nSVexpansion nSV nSVError C F1 MCC]
%         where
%             Err_Rate=error rate
%             AUC= area under the ROC curve
%             Err_RateA=error rate anomalous groups
%             Err_RateN=error rate normal groups
%             nSVexpansion=number of support vector expansion
%             nSV=number of support vectors
%             nSVError=number of support vector errors
%             C=parameter
%             F1=F1 score
%             MCC= Mathews statistic
% the models tested were: SMDDCPP=1, SMDDDA=2, SMDDDA with spherical
% normalization=3, OCSMM=4 and SVDD=5, the cell
% statistics(EPOCH, model, C, lambda), contains the statistics for each model, 
% for C =[1:1:10] corresponding {[1 1./(nTraining*[0.1:0.1:0.9])]} 
% and lambda=[1:1:4] corresponding {[sigma1 1/sigma1 1/sigma2 1/sigma3]}, 
% for example
% statistics(:,1,3,4), has the statistics for 200 runs, for model 1
% (SMDDCPP), C=3 (0.0110, this values makes 20% of the training set errors)
% and lambda = 1/sigma3,
% where sigma1=median(distances among all points in the dataset)
%       sigma2 = 0.1s quantile(distances among all points in the dataset)
%       sigma3 = 0.9 quantile(distances among all points in the dataset)


%run   /home/jorjasso/cvx/cvx_startup.m
run /home/jorjasso/Downloads/cvx/cvx_startup.m

experiementSloan(0)
experiementSloan(1)
experiementSloan(5)
experiementSloan(10)
experiementSloan(100)

%files 
%'experimentosSDSS0Sigma.mat'
%'experimentosSDSS1Sigma.mat'
%'experimentosSDSS2Sigma.mat'
%'experimentosSDSS3Sigma.mat'
%'experimentosSDSS4Sigma.mat'


%factor=3
%load (['experimentosSDSS' int2str(factor) 'Sigma.mat'])

%statistics=[Err_Rate AUC  Err_RateA Err_RateN nSVexpansion nSV nSVError C F1 MCC];

% matrixC is a matrix with statistics for some C =[1:1:10] correspond {[1 1./(nTraining*[0.1:0.1:0.9])]} 
% and lambda=[1:1:4] correspond {[sigma1 1/sigma1 1/sigma2 1/sigma3]}

C=1;lambda=1;
matrixC = cell2mat(statistics(:,:,C,lambda));


%plots
ind=mod(1:50,10);
figure
boxplot(100-matrixC(:,ind==1));           %ACC
figure
boxplot(matrixC(:,ind==2));               %AUC
figure
boxplot(100-matrixC(:,ind==3)); % ACC anomalous groups
figure
boxplot(100-matrixC(:,ind==4)); % ACC normal groups
figure
boxplot(matrixC(:,ind==5)); % number of support vector expansion
figure
boxplot(matrixC(:,ind==6)); % number of support vectors
figure
boxplot(matrixC(:,ind==7)); % number of support vectors errors
figure
boxplot(matrixC(:,ind==8)); % C
figure
boxplot(matrixC(:,ind==9)); % F1
figure
boxplot(matrixC(:,ind==0)); % MCC
figure
boxplot(lambdaValues); %lambda



%% PCA example, dimensionality reduction
load hald;
[coeff,score,latent,tsquare] = princomp(ingredients);
[coeff1,latent1,explained] = pcacov(cov(ingredients));

% is the zero mean centered data
data=bsxfun(@minus, ingredients,mean(ingredients));
[V,A]=eig(cov(data'*data));


% by theory data=score*coeff'
mean(mean(abs(data-score*coeff'))) % is small value
%score is the proyected data in PCA space, score=data*coeff
mean(mean(abs(score-data*coeff))) % is small value

%percent of the total variance for each PCA component
latent*100/sum(latent)

%dimensionality reduction in PCA dimension, taking 2 of 4 dimensions about
%97% of the variance
scoreRed=data*coeff(:,1:2)

%dimensionality reduction in input dimension, taking 2 of 4 dimensions
%97% of the variance
dataRed=score*coeff(1:2,:)'




%% image data
addpath ./anomalousScenes/SceneClass13/CALsuburb
addpath ./anomalousScenes/SceneClass13/MITcoast
addpath ./anomalousScenes/SceneClass13/MITforest
addpath ./anomalousScenes/Improved_Nystrom_Method

run   /home/jorjasso/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
addpath ./SVM-KM/






return


rmpath ./anomalousScenes/SceneClass13/


% extract features of all data 
    % for each image randomly compute 100patches, each patch, may be of
    % constant size (100x100 pixels) or randomly position and scale, 
    %(in reference, ), patch are determined as:
    % 1 evenly sampled grid: at 10x 10 pixels for image, the size of the
    % patch is randomly between scale of 10 to 30 pixels (this gives better results)
    % 2 randomly sampled: 500 randomly sampled patches for image, the size
    % is randomly sampled bewteen scale 10 to 30 pixels
    % 
    % each patch
    % gives a  128 dimensional SIFT feature vector
    % apply PCA to SIFT features to get 95 % of the variance
    % form the clusetr with those features

    % experiment    




    






%% conection with one class SVM OTRO PAPER, TRABAJAR EN ESTO PARA OTRO
%% PAPER
run /home/jorjasso/Downloads/cvx/cvx_startup.m
clc;clear;
load datosSameSigma
%experimentar con varios kappas y anchos de gaussianas para ver que pasa
%hacer un grafico de SMDDDCCPA con el SMDDCCP y el OCSVM para otro paper
kappa=1*ones(n,1);
kernelOp=3;
par=2^(1)

%test-----------------------
[R1,c1,cNorm1,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);


[rho,c,cNorm,alph,Io,I1,I2,I3]=SMDDCCPASOCSVM(S,muTrain,Sigma,kappa,C,kernelOp,par)

[K_trace,~]=kernelPM(S(Io),muTrain(Io),Sigma(Io),S(Io),muTrain(Io),Sigma(Io),kernelOp,par);
R=sqrt((K_trace(1)+cNorm-2*rho)/kappa(1));


[R1 R abs(R1-R)]
[c1' c' abs(c1'-c')]
[cNorm1 cNorm abs(cNorm1-cNorm)]
%cNorm-2*rho
[rho -0.5*(R1^2*kappa(1)-K_trace(1)-cNorm1)]
%---------------------------------


%SMDDCCP
[R,c,cNorm,alph1,Io,I1,I2,I3]=SMDDCCP(S,muTrain,Sigma,kappa,C,kernelOp,par);
plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,0)
saveas(gcf,['copiar/figSMDDCCPsameSigma'],'fig');
print(gcf, '-dpng', '-r0', ['copiar/figSMDDCCP.png']);
close(gcf);

% ypred=3 SVDDCPM as OCSVM
%ypred=3;
[R,c,cNorm,alph,Io,I1,I2,I3]=SMDDCCPASOCSVM(S,muTrain,Sigma,kappa,C,kernelOp,kernelParam)

rho=R;

[K_trace,~]=kernelPM(S(Io),muTrain(Io),Sigma(Io),S(Io),muTrain(Io),Sigma(Io),kernelOp,kernelParam);

R=sqrt((K_trace(1)+cNorm-2*rho)/kappa(1));

% plotContourKernelPM(S,muTrain,Sigma,R,cNorm,alph1,Io,I1,kernelOp,par,ypred)
% saveas(gcf,['copiar/figSMDDCCPAsOCSVM'],'fig');
% print(gcf, '-dpng', '-r0', ['copiar/figSMDDCCPAsOCSVM.png']);
% close(gcf);





%% DATASETS
%% Data with same covariance matrix
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSMDDDAvsSMDDDCPPModificados
n=10;
C=1;
L=100;
Sigma{2}=Sigma{2}*5.7;
muTrain(2,1)=-0.5
S{2} = mvnrnd(muTrain(2,:),Sigma{2},L);
muTrain(2,:)=mean(S{2});             % estimated mean
Sigma{2}=cov(S{2});

for i=1:n
    S{i}=mvnrnd(muTrain(i,:),Sigma{1},L);
    muTrain(i,:)=mean(S{i});             % estimated mean
    Sigma{i}=cov(S{1});
end

save datosSameSigma.mat


%% Data with same covariance matrix cov matrix is 0.5 times the datosSameSigma version
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSameSigma.mat
n=10;
C=1;
L=100;

for i=1:n
    S{i}=mvnrnd(muTrain(i,:),Sigma{1}*0.5,L);
    muTrain(i,:)=mean(S{i});             % estimated mean
    Sigma{i}=cov(S{1});
end

save datosSameSigmaII.mat

%% Data with same covariance matrix cov matrix is 0.25 times the datosSameSigma version
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSameSigma.mat
n=10;
C=1;
L=100;

for i=1:n
    S{i}=mvnrnd(muTrain(i,:),Sigma{1}*0.25,L);
    muTrain(i,:)=mean(S{i});             % estimated mean
    Sigma{i}=cov(S{1});
end

save datosSameSigmaIII.mat

%% Data with same covariance matrix cov matrix is 0.125 times the datosSameSigma version
run /home/jorjasso/Downloads/cvx/cvx_startup.m
load datosSameSigma.mat
n=10;
C=1;
L=100;

for i=1:n
    S{i}=mvnrnd(muTrain(i,:),Sigma{1}*0.125,L);
    muTrain(i,:)=mean(S{i});             % estimated mean
    Sigma{i}=cov(S{1});
end

save datosSameSigmaIV.mat

%% Sigma=U*V*U'=(U*sqrt(V))*(U*sqrt(V))'
N=100;
X=randn(N,2);
SS=cov(X);
%plot original data
figure; plot(X(:,1),X(:,2),'.r');
%scale and rotate data
scale=[3 0;0 1];
theta=30*pi/180;
rotate=[cos(theta) -sin(theta) ;sin(theta) cos(theta) ];
X=X*scale*rotate;
%plot
figure; plot(X(:,1),X(:,2),'.r');
grid on;
axis square equal;
% covariance
Sigma=cov(X);
%eigenvalue descomposition
[U V]=eig(Sigma);
%Sigma=U*V*U'=(U*sqrt(V))*(U*sqrt(V))'
for sd=[0.15  1  2]
    plot_gaussian_ellipsoid(mean(X), Sigma,sd,100,'k-')
end



close all
muTrain=[3 4; 7 9];
Sigma{1}=[4 0; 0 1];
Sigma{2}=[4 0; 0 1];
theta=30*pi/180;
R=real([cos(theta) -sin(theta) ;sin(theta) cos(theta) ]);
Sigma{2}=R*Sigma{2}*R';
for i=1:n
      %for sd=[1:0.7:3]
      %for sd=[0.15  0.5  1 ]
      for sd=[0.15    1 2]
      plot_gaussian_ellipsoid(muTrain(i,:), Sigma{i},sd,100,'k-')     
      end      
    %  [U D]=eig(Sigma{i});    
    % Plot first eigenvector
    %line([muTrain(i,1) muTrain(i,1)+sqrt(D(1,1))*U(1,1)],[muTrain(i,2) muTrain(i,2)+sqrt(D(1,1))*U(2,1)],'linewidth',3)
    % Plot second eigenvector
    %line([muTrain(i,1) muTrain(i,1)+sqrt(D(2,2))*U(1,2)],[muTrain(i,2) muTrain(i,2)+sqrt(D(2,2))*U(2,2)],'linewidth',3)
  end
  %%


%% training dataset
clear 
clc
for i=1:30
    n=poissrnd(100);
    S{i} = mvnrnd([1,1],eye(2),n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
%     figure
%     plot(S{i}(:,1),S{i}(:,2),'.r')
%     for sd=[0.15    1 2]
%       plot_gaussian_ellipsoid(mu(i,:), Sigma{i},sd,100,'k-')     
%     end      
end

%test set
scale=[1 0;0 1];
theta=45*pi/180;
rotate=[cos(theta) -sin(theta) ;sin(theta) cos(theta) ];
for i=1:10
    n=poissrnd(100);
    SVal{i} = S{i} *scale*rotate; muVal(i,:)=mean(SVal{i});SigmaVal{i}=cov(SVal{i});
%    figure
%     plot(SVal{i}(:,1),SVal{i}(:,2),'.g')
%     for sd=[0.15    1 2]
%       plot_gaussian_ellipsoid(muVal(i,:), SigmaVal{i},sd,100,'-g')     
%     end      
end

for i=11:20
    n=poissrnd(100);
    SVal{i} = mvnrnd([1,1],eye(2),n);muVal(i,:)=mean(SVal{i});SigmaVal{i}=cov(SVal{i});
%     figure
%     plot(S{i}(:,1),S{i}(:,2),'.r')
%     for sd=[0.15    1 2]
%       plot_gaussian_ellipsoid(mu(i,:), Sigma{i},sd,100,'k-')     
%     end      
end


%plot the means
figure;
plot(mu(:,1),mu(:,2),'.r')
hold on
plot(muVal(:,1),muVal(:,2),'.g')

training={S mu Sigma};
test={SVal muVal SigmaVal};

[sigma1 sigma2 sigma3]=findLambda(S)
kernelParam=1/sigma1
kappa=ones(length(S),1);
kernelOp=3;
ModelOp=1;
C=1;

yt=[ones(10,1);-ones(10,1)];
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
 
ModelOp=2;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
   0     1     0     0     4     4     0     1     1      1

ModelOp=3;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)

     0     1     0     0     3     3     0     1   1     1

ModelOp=4;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)

ModelOp=5;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)



%% 
%training
clear
clc
muGM = [-1.7 -1; 1.7 -1; 0 2];
I=0.2*eye(2);
sigmaGM = I;
p1=[0.33,0.33,0.33];%weigths for normal groups

obj = gmdistribution(muGM,sigmaGM,p1);
for i=1:50
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end

%test
scale=[1 0;1 0];
theta=45*pi/180;
rotate=[cos(theta) -sin(theta) ;sin(theta) cos(theta) ];
sigmaGM1=cov(S{1} *rotate);
sigmaGM2=cov(S{1} *scale);
sigmaGM3=cov(S{1} *scale*rotate);

obj = gmdistribution(muGM,cat(3,sigmaGM1,sigmaGM2,sigmaGM3),p1);
for i=1:20
    n=poissrnd(100);
    SVal{i} = random(obj,n); muVal(i,:)=mean(SVal{i});SigmaVal{i}=cov(SVal{i});
%   figure
%plot(SVal{i}(:,1),SVal{i}(:,2),'.g')
%     for sd=[0.15    1 2]
%       plot_gaussian_ellipsoid(muVal(i,:), SigmaVal{i},sd,100,'-g')     
%     end      
end

obj = gmdistribution(muGM,sigmaGM,p1);
for i=21:40
    n=poissrnd(100);
    SVal{i} = random(obj,n); muVal(i,:)=mean(SVal{i});SigmaVal{i}=cov(SVal{i});
%    figure
%     plot(SVal{i}(:,1),SVal{i}(:,2),'.g')
%     for sd=[0.15    1 2]
%       plot_gaussian_ellipsoid(muVal(i,:), SigmaVal{i},sd,100,'-g')     
%     end      
end

%plot the means
figure;
plot(mu(:,1),mu(:,2),'.r')
hold on
plot(muVal(:,1),muVal(:,2),'.g')

%
 [training,  test]=getDistributionBasedData(30,40,0);

S=training{1};
[sigma1 sigma2 sigma3]=findLambda()
kernelParam=1/sigma1
kappa=ones(length(S),1);
kernelOp=3;
ModelOp=1;
C=1;

yt=[ones(20,1);-ones(20,1)];
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
 
ModelOp=2;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
 

ModelOp=3;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)



ModelOp=4;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
20.0000    1.0000         0   40.0000    2.0000 2.0000    1.0000    1.0000    0.8333    0.6547

ModelOp=5;
statistics=prediction_MONQP(training, test, kernelParam,C,kappa,kernelOp,ModelOp,yt)
45.0000    0.5050   80.0000   10.0000    3.0000    3.0000         0    1.0000    0.3077    0.1400


