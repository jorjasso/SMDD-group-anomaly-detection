function [nonAnomalous, anomalousFirstType,anomalousSecondType,anomalousThirdType]=getPointBasedGroupAnomalies(N,OpPrint)
%  This function get: N non anomalous groups, N anomalousgroups (three different settings each)
% Point-based anomalies: agrregation of anomalous points

%%%%%%%%%%%%%%%%%
%NON-ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%

muGM = [-1.7 -1; 1.7 -1; 0 2];
I=0.2*eye(2);
sigmaGM = I;
p1=[0.33,0.64,0.03];%weigths for normal groups
p2=[0.33,0.03,0.64];%weigths for normal groups


obj = gmdistribution(muGM,sigmaGM,p1);
S=[];mu=[];Sigma=[];
for i=1:floor(N\2)
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
obj = gmdistribution(muGM,sigmaGM,p2);
for i=(floor(N\2)+1):N
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end

nonAnomalous={S mu Sigma};
%%%%%%%%%%%%%%%%%
%ANOMALOUS GROUPS
%%%%%%%%%%%%%%%%%

%First type of group anomalies
S=[];mu=[];Sigma=[];
for i=1:N
    n=poissrnd(100);
    S{i} = mvnrnd([-0.4,1],eye(2),n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end

anomalousFirstType={S mu Sigma};

%Second type of group anomalies
S=[];mu=[];Sigma=[];
p=[0.1 ,0.08,0.07  0.75];
obj = gmdistribution([muGM; 0.6 -1],sigmaGM,p);
for i=1:N
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
    %scatter(STest{i}(:,1),STest{i}(:,2),300,'.')
end

anomalousSecondType={S mu Sigma};

%Third type of group anomalies

p=[0.14  0.1 0.28 0.48];
obj = gmdistribution([muGM; -0.5 1],sigmaGM,p);
for i=1:N
    n=poissrnd(100);
    S{i} = random(obj,n);mu(i,:)=mean(S{i});Sigma{i}=cov(S{i});
end
anomalousThirdType={S mu Sigma};

if OpPrint==1
    for i=1:N
        set(gca,'FontSize',30,'fontWeight','bold')
        set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
        scatter(S{i}(:,1),S{i}(:,2),300,'.')
        %set the right path here
        saveas(gcf,['figuresGroupAnomalousDetectionMultimodal/figGMM' int2str(i)],'fig');
        print(gcf, '-dpng', '-r0', ['figuresGroupAnomalousDetectionMultimodal/figGMM' int2str(i) '.png']);
    end
end
