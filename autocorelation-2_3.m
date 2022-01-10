clear all


file_name=('signal.txt');
data=load(file_name);

n=data(:,1);
t=data(:,2);
U=data(:,3);
V=data(:,4);
W=data(:,5);


T = max(t);
tau = T/10;

%Bruit blanc

for i= 1:4

bb1 = wgn(1,length(t),-20*i);
bb2 = wgn(1,length(t),-20);

r=corrcoef(bb1,bb2)

[Rxx,lags]= xcorr(bb1,'normalized');


figure (1)
subplot(4,1,i);plot(lags,Rxx)
xlabel('tau')
ylabel('auto-corrélation')
title('dB=',i*-20)

end

figure (2)
subplot(2,1,1);plot(t,bb1)
title('Bruit Blanc')
xlabel('temps')
ylabel('Amplitude')
grid on

subplot(2,1,2);plot(lags,Rxx)
xlabel('tau')
ylabel('auto-corrélation')

grid on


% % Signal sinusoidal
% 
% Q= 0.1:0.01:0.9;
% 
% % for i=1:81
% 
% %tau(i)= Q(i) .* T./10;
% f = 0.01;          % fréquence
% A = 1;              % amplitude
% 
% C1= A*cos(2*pi.*t*f);
% C2= A*cos(2*pi.*(t+tau)*f);
% % noise=rand(size(t));
% % C=C+noise;
% 
% [Rxx,lags]= xcorr(C1,C2);
% 
% % r=corrcoef(C1,C2);
% 
% % o(i)=r(2,1);
% 
% 
% 
% figure (2)
% subplot(2,1,1);plot(t,C1)
% title('Signal sinusoidal')
% xlabel('temps')
% ylabel('Amplitude')
% grid on
% 
% subplot(2,1,2);plot(lags,Rxx)
% xlabel('tau')
% ylabel('auto-corrélation')
% % title('tau=',tau(i))
% grid on

% end

% figure (4)
% plot(tau,o,'-x')
% xlabel('tau')
% ylabel('AC coef')
% grid on

% Langevin
% 
% for y=1:5
% dt=0.0050;
% ll1= Langevin(0,1,tau,dt,length(t));
% ll2= Langevin(0,1,tau,dt,length(t));
% 
% [Rxx,lags]= xcorr(ll1,'unbiased');
% 
% r=corrcoef(ll1,ll2);
% 
% tau_intlg = trapz(lags,Rxx)
% tau_intlg = trapz(lags,Rxx)
% 
% figure (y)
% subplot(2,1,1);plot(t,ll1)
% title('Langevin')
% xlabel('temps')
% ylabel('Amplitude')
% grid on
% 
% 
% subplot(2,1,2);plot(lags,Rxx)
% xlabel('tau')
% ylabel('auto-corrélation')
% grid on
% end



function X = Langevin(Xmean, Xvar, T, dt, N)
	
	%return a signal given by the Langevin process
	% with:
	% * Xmean: the mean of the process
	% * Xvar: its variance
	% * T:its correlation time 
	% * dt: the time step
	% and N the number of time step
	
	dt_adim=dt/T;
	h=sqrt(Xvar*dt_adim);
	
	X=zeros(N,1);
	X(1)=randn()*sqrt(Xvar);
    for i=2:N
		dx = -(X(i-1) - Xmean) * dt_adim;
		dx = dx + randn()* h;
		X(i) = X(i-1) + dx ;
    end
	
end




