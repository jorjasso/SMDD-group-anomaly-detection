function [Err_Rate, Err_RateA,Err_RateN,ConfMat,AUC,stat]=Error_count(yt,predict_label)
% statistics for yt (true labels) vs predic_label
%output?  Err_Rate = error , ConfMat = confusion matrix,
% AUC = area under the curve ROC, stat =vector with the general estatistics:
% stat =[TruePositive TrueNegative FalsePositive FalseNegative TruePostiveRate SPeCifity PrecisionofPositiveValues NegativePredictiveValue FalsePositiveRate FalseDiscoveryRate ACCuracy F1score MatthewsCorrelationCoefficient];


%error for anomalous groups or class 1
tt=length(find(yt==1));
ytP=yt(1:tt); predict_labelP=predict_label(1:tt);
ConfMat=zeros(2,2);
ConfMat(1,2)=length(find((ytP==1)&(sign(predict_labelP)==-1)));
ConfMat(2,1)=length(find((ytP==-1)&(sign(predict_labelP)==1)));

nt=length(ytP);
Err_RateA=100*(ConfMat(1,2)+ConfMat(2,1))/nt;

%error for normal groups or class -1

ttN=length(find(yt==-1));
ytN=yt(tt+1:tt+ttN); predict_labelN=predict_label(tt+1:tt+ttN);
ConfMat=zeros(2,2);
ConfMat(1,2)=length(find((ytN==1)&(sign(predict_labelN)==-1)));
ConfMat(2,1)=length(find((ytN==-1)&(sign(predict_labelN)==1)));

nt=length(ytN);
Err_RateN=100*(ConfMat(1,2)+ConfMat(2,1))/nt;

%ROC ANALISYS
testpred=predict_label;
testt=yt;
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
%ROC curve
%figure(1);hold off
cspec = 1-spec; %FPR
cspec = cspec(end:-1:1);
sens = sens(end:-1:1);
%plot(cspec,sens,'k')
%title('ROC curve')
AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));

%[X,Y,T,AUC] = perfcurve(yt,predict_label,1)
%plot ROC
%plot(X,Y)
%xlabel('False positive rate'); ylabel('True positive rate')
%title('ROC for classification by logistic regression')

%general statistics

%confusion matrix
ConfMat=zeros(2,2);
ConfMat(1,1)=length(find((yt==1)&(sign(predict_label)==1)));
ConfMat(1,2)=length(find((yt==1)&(sign(predict_label)==-1)));
ConfMat(2,1)=length(find((yt==-1)&(sign(predict_label)==1)));
ConfMat(2,2)=length(find((yt==-1)&(sign(predict_label)==-1)));

%error Total
nt=length(yt);
Err_Rate=100*(ConfMat(1,2)+ConfMat(2,1))/nt;

TP=ConfMat(1,1);%true positive: normal grous correctly detected as normal
TN=ConfMat(2,2);%true negative: outlier detected as outlier look this for number of outliers detected
FN=ConfMat(1,2);%false negative: normal groups detected as outlier, TYPE II ERROR (table for wikipedia is rotated 90 degrees :))
FP=ConfMat(2,1);%false positive: outliers detected as normal,  TYPE I ERROR
TPR=TP/(TP+FN);%sensitivity or true positive rate or recall,    % useful for ROC
SPC=TN/(FP+TN);%specivicity or true negative rate, look this for outlier detection, because we are looking for the negatives (outliers)?
PPV=TP/(TP+FP);%presicion or positive predictive value
NPV=TN/(TN+FN);%negative predictive value
FPR=FP/(FP+TN);%fall-out or false positive rate or 1-SPC  % useful for ROC
FDR=FP/(TP+FP);%false discovery rate or 1-PPV
ACC=(TP+TN)/(TP+FN+FP+TN);%accuracy
F1=2*TP/(2*TP+FP+FN);%F1 score
MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));%Matthews correlation coefficient

stat=[TP TN FP FN TPR SPC PPV NPV FPR FDR ACC F1 MCC];