% USVDD prediction without kernels
% input
%            Xtrain,trSigmaTrain,kappaTrain = uncertain training set
% output
%            y= predictions {-1,1} belong or not belong to the support

function y=predictUSVDD(XTest,trSigmaTest,kappaTest,Xtrain,trSigmaTrain,kappaTrain,C,op)

switch lower(op)
    case {'one'}
        disp('Metric one')
        %training
        [R,c,alph1,Io,Isurface,Ierrors,IinSphere]=USVDD(Xtrain,trSigmaTrain,kappaTrain,C);
        %prediction
        A=bsxfun(@minus, XTest, c); % each row is x_i-c        
        y=sign(bsxfun(@minus, R^2, diag(A*A'))); %    each row is ||x_i-c||^2
    case 'two'
        disp('Metric two')
    case 'three'
        disp('Metric three')
    otherwise
        disp('Other')
end%