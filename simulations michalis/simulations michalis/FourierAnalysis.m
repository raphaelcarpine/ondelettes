%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FOURIER ANALYSIS:
%
% created by Michalis Fragiadakis, Dec 2013
% mfrag@mail.ntua.gr
%
% please report any bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
includeFile

npts = acc.nsteps;
asignal = acc.rec;
dt = acc.dtacc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOURIER ANALYSIS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=npts;
FASsignal=fft(asignal,N)*dt;
dw=2*pi/(N*dt);
w=linspace(0,(N-1)*dw,N);
f=w/(2*pi);
nyq=ceil(N/2);

% plot the FAS
figure()
hold on; grid on; box on;
plot(w(1:nyq),abs(FASsignal(1:nyq)),'r-')
xlabel('frequency (Hz)','FontSize',16)
ylabel('Fourier Amplitude','FontSize',16)
set(gca,'xscale','log','xlim',[0,w(nyq)],'FontSize',14)

% direct calculation - without the fft command of Matlab
n=1:npts;
for k=1:nyq
    fo(k)=dt*sum(asignal'.*exp( -i*(k*dw).*(n*dt) ));
end;
% plot((1:nyq)*dw,abs(fo(1:nyq)),'b--')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INVERSE FOURIER ANALYSIS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% get back the signal combing FAS with the phase

fi=angle(FASsignal(nyq:N-1));
%fi   = angle(FASsignal(1:nyq));
z    = FASsignal(1:nyq).*exp(i*fi);
atmp = (1/dt)*ifft(z,N,'symmetric');

% plot the new record and compare with the initial
figure()
hold on; box on; grid on;
plot(dt*(1:npts),-asignal,'b-','LineWidth',2)
plot(dt*(1:npts),atmp,'r-','LineWidth',1)
xlabel('time (sec)','FontSize',20)
ylabel('acceleration (g)','FontSize',20)
hleg = legend('initial signal','inverse FFT');
set(hleg,'FontSize',20)


% direct calculation - without the ifft command of Matlab
k=1:nyq;
for n=1:npts
    tmp = fo(1:nyq) .* exp(i.*(k*dw)*(n*dt)); %<- warning, use either i or -i
    aa(n)= (1/pi)*dw * sum(tmp);
end;
% plot((1:npts)*dt,-aa,'g--')