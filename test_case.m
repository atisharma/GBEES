% MATLAB file
% loads and plots P

function [P,X,Y,V,V0,Veps,Vfft, Vhat]=test_case(dx);
	%addpath('./export_fig');
	
	P0=load('-ascii','0.00_LeVeque_test_case.asc');
	Peps=load('-ascii','6.28_final_eps_LeVeque_test_case.asc');
	P=load('-ascii','6.28_final_pure_LeVeque_test_case.asc');

	L=1; dy=dx;
	X=[-L:dx:L]; Y=X; N=length(X);
	[XG,YG] = meshgrid(X,Y);
	
	V0 = griddata(P0(:,1),P0(:,2),P0(:,7),XG,YG); 
	Veps = griddata(Peps(:,1),Peps(:,2),Peps(:,7),XG,YG); 
	V = griddata(P(:,1),P(:,2),P(:,7),XG,YG); 

	% analytical diffusion
	mu=4e-5; % based on minimising KL to test case dx=.01, dt=.001
	Vfft = V0;
	Vfft(find(isnan(Vfft)))=0; % set not defined p to zero
	M=size(Vfft,1);
	N=size(Vfft,2);
	Vhat = fft2(Vfft);
	
	% convert DFT wavenumber to actual wavenumbers
	% Thank you Krogstad
	kx1 = mod( 1/2 + (0:(M-1))/M , 1 ) - 1/2;
	kx = kx1*(2*pi/dx);
	ky1 = mod( 1/2 + (0:(N-1))/N , 1 ) - 1/2;
	ky = ky1*(2*pi/dy);
	[KX,KY] = meshgrid(kx,ky);
	
	D = exp(-mu *2*pi* (KX.*KX + KY.*KY)); % differential operator in Fourier space
	Vhat = D.*Vhat;
	Vfft = real(ifft2(Vhat,'symmetric'));
	
	disp(['mu=' num2str(mu)])
	disp(['D_KL(P0, [P_eps P P_diffusion] )=' ...
		num2str([D_KL(V0,Veps) D_KL(V0,V) D_KL(V0,Vfft)]) ]);
	disp(['D_KL(P_eps, P_diffusion )=' num2str([D_KL(Veps,Vfft)]) ]);
	disp(['D_KL(P, P_eps )=' num2str([D_KL(Veps,V)]) ]);
	

	figure(1);clf;
	subplot(1,4,1)
	surf(Y,X,V0,'linestyle','none');
	axis equal

	subplot(1,4,2)
	surf(Y,X,Veps,'linestyle','none');
	axis equal

	subplot(1,4,3)
	surf(Y,X,V,'linestyle','none');
	axis equal

	subplot(1,4,4)
	surf(X,Y,Vfft,'linestyle','none');
	axis equal
	
	figure(4);clf
	ix=find(X==sqrt(min(X.*X))); % x=0 plane
	plot(X,Veps(ix,:),'b--','linewidth',1.5); % truncation 1e-3
	hold on
	%plot(X,V(ix,:),'c--'); % almost no truncation
	plot(X,V0(ix,:),'k','linewidth',1.5); % analytic
	plot(X,real(Vfft(ix,:)),'r-.','linewidth',1.5); % analytic + diffusion
	ylim([0 1.1])
	xlim([-1 1])
	xlabel('y')
	%axis equal
	%axis tight
	set(gcf, 'PaperPosition', [0 0 6 3]);
	print -depsc2 x0_new.eps

return;
