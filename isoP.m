% MATLAB file
% loads and plots P
% pass the dx you want to use in the plot
% Eg if you used dx=0.4 in the GBEES code, try 0.8 here (for faster plot)
% or 0.4 for more accurate

function []=isoP(dx);
	addpath('./export_fig');
	edge1=0.005;
	edge2=0.0005;
	edge3=0.00005;


%%%%%%%%%% Time-stepped plot
	figure(1);clf;
	truth=load('-ascii','build/Release/truth.asc');
	frame=load('-ascii','build/Release/frame.asc');
	L=1.2*max([max(abs(frame(:,1))) max(abs(frame(:,2))) max(abs(frame(:,3)))]);
	X=[-L:dx:L]; Y=X;Z=X; N=length(X);
	
	c1=[edge1 (1-edge1) 1]; c2=	[edge2 (1-edge2) 1]; c3= [edge3 (1-edge3) 1];
	
	for t=0:.2:1
		eval(['plot_distrib(pload(''build/Release/' num2str(t,'%0.2f') '_distribution.asc''),X,Y,Z,c1,c2,c3);']);
	end
	
	plot3(frame(:,1),frame(:,2),frame(:,3),'g-', 'linewidth',.5,'linesmoothing','on');
	lighting phong; drawnow; view(-109,14);  hold on;
	light('Position',[-1 0 0]);
	plot3(truth(:,1),truth(:,2),truth(:,3),'k-','linesmoothing','on','linewidth',2);
	plot3(truth(end,1),truth(end,2),truth(end,3),'k*');

	view(-109,14)
	axis([-L L -L L -L*2 L]);
	axis equal

	set(gcf,'color','none')
	set(gca,'color','none')
	export_fig(['Lorenz_GBEES_transparent.png'],'-png');
	set(gcf,'color',[1 1 1])
	export_fig(['Lorenz_GBEES.png'],'-png');
	view(-109,14); print('-depsc2','-opengl','-r600',['Lorenz_GBEES_1.eps']);
	view(-31,2);   print('-depsc2','-opengl','-r600',['Lorenz_GBEES_2.eps']);
	drawnow

return

%%%%%%%%%% Time-stepped plot with swarm
	figure(2);clf;
	
	for t=0:.2:1
		eval(['plot_distrib(pload(''build/Release/' num2str(t,'%0.2f') '_distribution.asc''),X,Y,Z,c1,c2,c3);']);
	end

	for n=1:200
		bee=load('-ascii',['build/Release/swarm_' num2str(n) '.asc']);
		plot3(bee(:,1),bee(:,2),bee(:,3),'g-', 'linewidth',.5,'linesmoothing','on');hold on
		for t=0:.2:1
			ti=find(abs(bee(:,4)-t)==min(abs(bee(:,4)-t)));
			plot3(bee(ti,1),bee(ti,2),bee(ti,3),'k+');
		end
	end
	plot3(truth(:,1),truth(:,2),truth(:,3),'k-','linesmoothing','on','linewidth',2);
	plot3(truth(end,1),truth(end,2),truth(end,3),'k*');
	
	lighting phong; drawnow; view(-109,14);
	light('Position',[-1 0 0]);
	view(-109,14)
	axis([-L L -L L -L*2 L]);
	axis equal

	set(gcf,'color','none')
	set(gca,'color','none')
	export_fig(['Lorenz_GBEES_swarm_transparent.png'],'-png');
	set(gcf,'color',[1 1 1])
	export_fig(['Lorenz_GBEES_swarm.png'],'-png');
	view(-109,14); print('-depsc2','-opengl','-r600',['Lorenz_GBEES_swarm_1.eps']);
	view(-31,2);   print('-depsc2','-opengl','-r600',['Lorenz_GBEES_swarm_2.eps']);
	drawnow

return

%%%%%%%%%% Pre & post measurement plot
	figure(3);clf;
	plot_distrib(pload('build/Release/1.00_distribution_pre_meas.asc'),X,Y,Z,c1,c2,c3);
	c1=[1 edge1 (1-edge1)]; c2=	[1 edge2 (1-edge2)]; c3= [1 edge3 (1-edge3)];
	plot_distrib(pload('build/Release/1.00_distribution_post_meas.asc'),X,Y,Z,c1,c2,c3);
	
	plot3(frame(:,1),frame(:,2),frame(:,3),'g-', 'linewidth',.5,'linesmoothing','on');
	lighting phong; drawnow; view(-109,14);
	light('Position',[-1 0 0]);
	plot3(truth(:,1),truth(:,2),truth(:,3),'k-','linesmoothing','on','linewidth',2);
	plot3(truth(end,1),truth(end,2),truth(end,3),'k*');

	view(-109,14)
	axis([-L L -L L -L*2 L]);
	axis equal

	set(gcf,'color','none')
	set(gca,'color','none')
	export_fig(['Lorenz_GBEES_transparent_prepost.png'],'-png');
	set(gcf,'color',[1 1 1])
	export_fig(['Lorenz_GBEES_prepost.png'],'-png');
	view(-109,14); print('-depsc2','-opengl','-r600',['Lorenz_GBEES_prepost_1.eps']);
	view(-31,2);   print('-depsc2','-opengl','-r600',['Lorenz_GBEES_prepost_2.eps']);
	drawnow

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function []=rotateP();
	for th=70:2:430
		view([th 30])
		drawnow
	end
return;

function P=pload(pstr);
	try
		P=load('-ascii',pstr);
	catch
		P=nan;
		warning(['Loading error; ' pstr]);
	end
return

function plot_distrib(P,X,Y,Z,c1,c2,c3);
%	try
		N=length(X);
		V=nan*single(zeros(N,N,N));
		offset=find(abs(X)==min(abs(X)));
		for n=1:length(P)
			j=find(abs(X-P(n,1))==min(abs(X-P(n,1))));
			i=find(abs(Y-P(n,2))==min(abs(Y-P(n,2))));
			k=find(abs(Z-P(n,3))==min(abs(Z-P(n,3))));
			V(i,j,k) = P(n,10);
		end

		% fixed isovalues
		edge1=0.005;
		edge2=0.0005;
		edge3=0.00005;
		
		p=patch(isosurface(X,Y,Z, V, edge1)); hold on
		set(p, 'EdgeColor', 'none','facecolor',c1,'facealpha',0.5,'linesmoothing','on');
		p=patch(isosurface(X,Y,Z, V, edge2));
		set(p, 'EdgeColor', 'none','facecolor',c2,'facealpha',0.5,'linesmoothing','on');
		p=patch(isosurface(X,Y,Z, V, edge3));
		set(p, 'EdgeColor', 'none','facecolor',c3,'facealpha',0.5,'linesmoothing','on');
%	catch
%		warning('Plotting error');
%	end
return;