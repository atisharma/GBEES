% MATLAB file
% loads and plots P

function []=surfP(dx);
	addpath('./export_fig');
	P=load('-ascii','build/Release/distribution.asc');
	truth=load('-ascii','build/Release/truth.asc');
	frame=load('-ascii','build/Release/frame.asc');


	L=3*max([max(abs(frame(:,1))) max(abs(frame(:,2)))]);

	X=[-L:dx:L]; Y=X; N=length(X);

	V=single(zeros(N,N));
	for n=1:length(P)
		j=max(find(X<P(n,1)));
		i=max(find(Y<P(n,2)));
		V(i,j) = P(n,7)+V(i,j);
	end

	for n=1:1
		figure(n);clf;
		contour(Y,X,V);hold on;
		plot(frame(:,1),frame(:,2),'color',[.5 .5 .4], 'linewidth',.001);
		plot(truth(:,1),truth(:,2),'r-','linesmoothing','on','linewidth',2);
		plot(truth(end,1),truth(end,2),'g*');
		set(gca,'visible','off')
		axis([-L L -L L]); axis equal
		set(gcf,'color','none')
%		export_fig(['2D_DBEES' num2str(n) '_transparent.png'],'-png');
		set(gcf,'color',[0 0 0])
%		export_fig(['2D_DBEES' num2str(n) '.png'],'-png');
	end

return;
