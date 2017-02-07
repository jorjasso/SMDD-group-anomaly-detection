function [nonAnomalous, anomalous]=getDistributionBasedData(N,OpPrint)
% This function get: N non anomalous groups, N anomalousgroups 
% Distribution-based anomalies

%%%%%%%%%%%%%%%%%
%NON-ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%

muGM = [-1.7 -1; 1.7 -1; 0 2];
I=0.2*eye(2);
sigmaGM = I;
p1=1/3*ones(3,1);%weigths for normal groups

obj = gmdistribution(muGM,sigmaGM,p1);
S=[];mu=[];Sigma=[];
for i=1:N
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
%     figure
%     plot(S{i}(:,1),S{i}(:,2),'.g')
%     pause
end

nonAnomalous={S mu Sigma};

%%%%%%%%%%%%%%%%%
%ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%

%First type of group anomalies

scale=[1 0;0 1];
theta=45*pi/180;
rotate=[cos(theta) -sin(theta) ;sin(theta) cos(theta) ];
sigmaGM1=cov(S{1} *rotate);
sigmaGM2=cov(S{1} *scale);
sigmaGM3=cov(S{1} *scale*rotate);
% group anomalies
obj = gmdistribution(muGM,cat(3,sigmaGM1,sigmaGM2,sigmaGM3),p1);
S=[];mu=[];Sigma=[];
for i=1:N
    n=poissrnd(100);
    S{i} = random(obj,n); mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
    %  figure
    %plot(STest{i}(:,1),STest{i}(:,2),'.g')
    %pause
    %     for sd=[0.15    1 2]
    %       plot_gaussian_ellipsoid(muTest(i,:), SigmaTest{i},sd,100,'-g')
    %     end
end

%plot the means
% figure;
% plot(mu(:,1),mu(:,2),'.r')
% hold on
% plot(muTest(:,1),muTest(:,2),'.g')

anomalous={S mu Sigma};

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
