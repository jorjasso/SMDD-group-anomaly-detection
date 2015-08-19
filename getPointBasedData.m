function [training, test]=getPointBasedData(nTraining,nTest,OpPrint)
%  References: Support measure data description support measure description for group anomaly detection
% jorge.jorjasso@gmail.com
% function [training, validation,test]=generateAnomalusGroupDetectionData(nTraining,nValidation,nTest)
% generate the training set, the validation set and the test set for the
% experiment, nTraining=number of groups for training,nValidation=number of
% groups for validation,nTest=number of groups for test
%OpPrint = 0, no print figure
%OpPrint = 1, print figure
muGM = [-1.7 -1; 1.7 -1; 0 2];
I=0.2*eye(2);
sigmaGM = I;
p1=[0.33,0.64,0.03];%weigths for normal groups
p2=[0.33,0.03,0.64];%weigths for normal groups


%50 no anomalous groups  with number of points per group N~poisson(100)
nFirstGroup=(48/100*nTraining); % 48% of the data

obj = gmdistribution(muGM,sigmaGM,p1);
for i=1:nFirstGroup
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
nSecondGroup=(52/100*nTraining); % 52% of the data
obj = gmdistribution(muGM,sigmaGM,p2);
for i=(nFirstGroup+1):(nFirstGroup+nSecondGroup)
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end



%Test data set
%-------------
STest=[];muTest=[];SigmaTest=[];
n1=ceil(nTest*33/100);
for i=1:n1
    n=poissrnd(100);
    STest{i} = mvnrnd([-0.4,1],eye(2),n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
end
%Ten group anomaly for point individually normal
n2=ceil(nTest*(33*0.5)/100);
p=[0.1 ,0.08,0.07  0.75];%weigths for anomalous groups
obj = gmdistribution([muGM; 0.6 -1],sigmaGM,p);
for i=(n1+1):(n1+n2)
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
    %scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
end
% 33*2 %of data
n3=ceil(nTest*(33*2)/100);
p=[0.14  0.1 0.28 0.48];%weigths for anomalous groups
obj = gmdistribution([muGM; -0.5 1],sigmaGM,p);
for i=(n1+n2+1):n3
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
end
%10 normal groups

obj = gmdistribution(muGM,sigmaGM,p1);
for i=(n3+1):(n3+5)
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
end
obj = gmdistribution(muGM,sigmaGM,p2);
for i=(n3+6):nTest
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
end


training={S mu Sigma};
test={STest muTest SigmaTest};

if OpPrint==1
    for i=1:nTest
        set(gca,'FontSize',30,'fontWeight','bold')
        set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
        scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
        saveas(gcf,['figuresGroupAnomalousDetectionMultimodal/figGMM' int2str(i)],'fig');
        print(gcf, '-dpng', '-r0', ['figuresGroupAnomalousDetectionMultimodal/figGMM' int2str(i) '.png']);
    end
end
