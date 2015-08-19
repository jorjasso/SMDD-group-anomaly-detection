function [training,  test]=generateAnomalusGroupDetectionData(nTraining,nTest,OpPrint)
% function [training, validation,test]=generateAnomalusGroupDetectionData(nTraining,nValidation,nTest)
% generate the training set, the validation set and the test set for the
% experiment, nTraining=number of groups for training,nValidation=number of
% groups for validation,nTest=number of groups for test
%OpPrint = 0, no print figure
%OpPrint = 1, print figure
muGM = [-1.7 -1; 1.7 -1; 0 2];
I=0.2*eye(2);
sigmaGM = I;
p=ones(1,3)/3;%weigths for normal groups
obj = gmdistribution(muGM,sigmaGM,p);

%Training data
%---------------
%50 no anomalous groups  with number of points per group N~poisson(100)
for i=1:nTraining
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end


%Test data set
%-------------
STest=[];muTest=[];SigmaTest=[];

%Ten group anomaly for point individually normal

p=[0.85,0.08,0.07];%weigths for anomalous groups
obj = gmdistribution(muGM,sigmaGM,p);
%for i=11:15

for i=1:10
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
    %scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
end
% 33*2 %of data

p=[0.04,0.48,0.48];%weigths for anomalous groups
obj = gmdistribution(muGM,sigmaGM,p);
%for i=16:20
for i=11:20
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
end
%30 normal groups

p=ones(1,3)/3;%weigths for normal groups
obj = gmdistribution(muGM,sigmaGM,p);
%for i=21:30
for i=21:nTest
    n=poissrnd(100);
    STest{i} = random(obj,n);muTest(i,:)=mean(STest{i}); SigmaTest{i}=cov(STest{i});
end

training={S mu Sigma};
test={STest muTest SigmaTest};

if OpPrint==1
    for i=1:nTest
        set(gca,'FontSize',30,'fontWeight','bold')
        set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
        scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
        saveas(gcf,['figuresGroupAnomalousDetection/figGMM' int2str(i)],'fig');
        print(gcf, '-dpng', '-r0', ['figuresGroupAnomalousDetection/figGMM' int2str(i) '.png']);
    end
end