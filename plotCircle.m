% PLOTS ellipsoildal data and, USVDD solution for those data
% return the handle for the center and the radio plots
% opt=0 no labels in the data, opt=1 for labels in the data
% colorAndType=color y tipo de linea
% lineWidth=line width
% Example:  [h1,h11]=plotCircle(X,Sigma,R,c,0,'--r',10)
function [h1,h11]=plotCircle(X,Sigma,R,c,opt,colorAndType,lineWidth)

addpath ./exportFig

[n,d]=size(X)
plot(X(:,1),X(:,2),'or')
for i=1:n
    for sd=[1:1:2]
    plot_gaussian_ellipsoid(X(i,:), Sigma{i},sd,100,'b-',gca)    
    end
    %plot the axis given by the eigen vectors
    %[U D]=eig(Sigma{i});    
    % Plot first eigenvector
    %line([X(i,1) X(i,1)+sqrt(D(1,1))*U(1,1)],[X(i,2) X(i,2)+sqrt(D(1,1))*U(2,1)],'linewidth',3)
    % Plot second eigenvector
    %line([X(i,1) X(i,1)+sqrt(D(2,2))*U(1,2)],[X(i,2) X(i,2)+sqrt(D(2,2))*U(2,2)],'linewidth',3)
end

hold on
%plotting circle uncertain SVDD
xc=c(1);yc=c(2)
ang=0:0.01:2*pi;
x=R*cos(ang)+xc;
y=R*sin(ang)+yc;
h1=plot(x,y,colorAndType,'LineWidth',lineWidth);
if opt==1
    labels = cellstr( num2str([1:n]') );  %' # labels correspond to their order
    text(X(:,1), X(:,2), labels, 'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right')
end
h11=plot(xc,yc,'xb')


%xlabel('X(1)','FontSize',14)
%ylabel('X(2)','FontSize',14)

%legend
%h_legend = legend([h1 h11 ],nameMethod,nameCenter);
%set(h_legend,'FontSize',14);
%title('\it{USVDD and SVDD solutions}','FontSize',14)


%%%%% produce high quality figures with export_fig
%%%%% https://sites.google.com/site/oliverwoodford/software/export_fig

 


