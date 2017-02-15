function th=getThreshold(alph1)
%function th=getThreshold(alph), th is the threshold to compute
% 0<alpha<C due to numerical precision by th<alph <C
%
th=eps;

% nSV=[];
% iz=length(alph1);
% der=0;
% bani=0;
% band=0;
% j=1;
% thresholds=2.^[-100:1:1];
% %computing left and righ limits [iz der] for the nSVues 
% for i=thresholds
%     nSV=[nSV,length(find(alph1>i))];
%     if iz>nSV(j)&bani==0
%         iz=j;
%         bani=1;
%     end
%     if der==nSV(j)&band==0
%         der=j;
%         band=1;
%     end
%     j=j+1;
% end
% 
% 
% for j=iz:der
%     if j+1>der
%         if  length(nSV(j))<length(nSV(j+1))
%             th = thresholds(j);
%         else
%             th=0.000009;
%         end
%     end
% end


 



