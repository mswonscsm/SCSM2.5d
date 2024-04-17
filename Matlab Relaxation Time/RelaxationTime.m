% RelaxationTime script calculates the relaxation time of stress and strain
% for input file: relaxation_time.inp
% Methodlogy is based on Emmerich and Korn, Incorporation of attenuation
% into time-domain computations of seismic wave fields, Geophysics, 1987,
% 52:9, 1252-1264.
% Strain Relaxation Time: tau_eta_e
% Stress Relaxation Time: tau_sig_e

clc; clear all
freq=input('Central Frequency ?');
f_ratio=101;
fmin0=1;
fmin=exp(log(freq)-log(12)/2);
fmax=12*fmin;
wmin=2*pi*fmin;%/sqrt(f_ratio);
wmax=2*pi*fmax;%*sqrt(f_ratio); 
xfrq=[wmin:(wmax-wmin)/(10000-1):wmax]; %w
x_axis=xfrq/(2*pi);
t=[0:10^-4:0.5];
Legend=cell(1,1); ii=0;
%%
%Constant Q
QPval=zeros(1,length(xfrq));
Q=input('Attenuation factor Q ?');
QPval(:)=Q;
%%
ns=input('How many springs ?');
nsn=size(ns,2);
y2=zeros(size(ns,2),size(t,2));
ii=0;
for N=ns
    ii=ii+1;
    Legend{ii}=strcat(num2str(N));

        [tau_eta_e, tau_sig_e]=Emmerich(f_ratio,QPval,wmin,wmax, xfrq, N);
        ex_e=zeros(N,size(t,2));
        ex_e=exp(-t./tau_sig_e');
        ex_e=((1-tau_eta_e./tau_sig_e))./N*ex_e; %N ->sum
        x=t;
        y=2.25*10^9*(1-ex_e);
        y2(ii,:)=y;

    figure(2),hold on,plot(t,y)
    pbaspect([2 1 1])
end

ylabel('Bulk Modulus (GPa)')
xlabel('t(s)')

y_log=log(y2);
theta_e=1./tau_sig_e;
kl_e=(tau_eta_e./tau_sig_e-1)/N;
%%
Qval=zeros(1,length(xfrq));
Qval(:)=QPval;
save('Qval.mat','xfrq','Qval')
x=[sqrt(abs(theta_e)) sqrt(abs(kl_e))];

%%
Q_f1=QpDraw(QPval, xfrq, tau_eta_e,tau_sig_e);

%%
figure(4),semilogx(xfrq/(2*pi),1./Qval, 'k')
hold on
figure(4),semilogx(xfrq/(2*pi),Q_f1, '--')
hold on

xlabel('log10(f)'); ylabel('1/Q')
legend('1/Qref','SolvOpt')

pbaspect([2 1 1])
