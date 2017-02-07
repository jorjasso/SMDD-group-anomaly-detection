function [nonAnomalosDataset, anomalousDataset]=getDatasets(nTraining,nTest,OpPrint)
%   input
%        type
%        N= number of observations for a non-anomalous dataset
%        NA = number of observations for a anomalous dataset
%        OpPrint = {0,1}, no print or print a figure
%   output
%        [non-anomalous dataset, anomalous dataset]


% non-anomalous dataset (Xiong setup)
%-----------------------
muGM = [-1.7 -1; 1.7 -1; 0 2];
sigmaGM=0.2*eye(2);
p1=1/3*ones(3,1);%weigths for normal groups
p2=[0.84,0.08,0.08]; p2(3)=1-p2(1)-p2(2);

obj = gmdistribution(muGM,sigmaGM,p1);
for i=1:nTraining
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
%     figure
%     plot(S{i}(:,1),S{i}(:,2),'.g')
%     pause
end

nonAnomalosDataset={S mu Sigma}

% anomalous dataset
%-----------------------

%scale=[1 0;1 0];
scale=[1 0;0 1];
theta=45*pi/180;
rotate=[cos(theta) -sin(theta) ;sin(theta) cos(theta) ];
sigmaGM1=cov(S{1} *rotate);
sigmaGM2=cov(S{1} *scale);
sigmaGM3=cov(S{1} *scale*rotate);
% group anomalies
obj = gmdistribution(muGM,cat(3,sigmaGM1,sigmaGM2,sigmaGM3),p1);
for i=1:floor(nTest/2)
    n=poissrnd(100);
    STest{i} = random(obj,n); muTest(i,:)=mean(STest{i});SigmaTest{i}=cov(STest{i});
    %  figure
    %plot(STest{i}(:,1),STest{i}(:,2),'.g')
    %pause
    %     for sd=[0.15    1 2]
    %       plot_gaussian_ellipsoid(muTest(i,:), SigmaTest{i},sd,100,'-g')
    %     end
end

anomalousDataset={STest muTest SigmaTest}

%plot the means
% figure;
% plot(mu(:,1),mu(:,2),'.r')
% hold on
% plot(muTest(:,1),muTest(:,2),'.g')

if OpPrint==1
    for i=1:nTest
        figure
        set(gca,'FontSize',30,'fontWeight','bold')
        set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
        scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
        
        saveas(gcf,['figuresGroupAnomalousDetection/figGMM' int2str(i)],'fig');
        print(gcf, '-dpng', '-r0', ['figuresGroupAnomalousDetection/figGMM' int2str(i) '.png']);
        
        
    end
end
