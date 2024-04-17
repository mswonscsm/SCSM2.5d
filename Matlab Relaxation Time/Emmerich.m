function [tau_eta, tau_sig]=Emmerich(f_ratio,QPval,wmin,wmax, xfrq, N)
K=2*N-1;

theta = wmin*(wmax/wmin).^([0:N-1]/(N-1));
wk=wmin*(wmax/wmin).^([0:K-1]/(K-1));
QPval_K=interp1(xfrq, QPval, wk);
over_Q=1./QPval_K;
%%
AK=zeros(K,N);
QB=zeros(K,1);
QB(:)=over_Q;
for k=1:K
    for l=1:N
        AK(k,l)=(wk(k)*(theta(l)-wk(k)*QB(k)))/(theta(l)^2+wk(k)^2);
    end
end
kl=AK\QB;
%%
Q_freq1=zeros(1,size(xfrq,2));
Q_freq2=zeros(1,size(xfrq,2));
for l=1:N
    Q_freq1=Q_freq1+xfrq.*theta(l)*kl(l)./(theta(l)^2+xfrq.^2);
    Q_freq2=Q_freq2+xfrq.^2*kl(l)./(theta(l)^2+xfrq.^2);
end
Q_f=Q_freq1./(1+Q_freq2);

%figure(3),hold on,plot(log10(xfrq/(2*pi)),ones(length(xfrq),1)./QPval)
%figure(3),semilogx(wk./(2*pi),QB, xfrq/(2*pi),Q_f, '--') 
%figure(3),hold on,semilogx(log10(xfrq/(2*pi)),Q_f,'b')
%legend('Q-ref', 'Yang', 'Emmerich')
%xlabel('log10(f)'); ylabel('1/Q')

%er=sum(1./QPval-Q_f).^2
%%
tau_sig=1./theta;
tau_eta=(1+N.*kl')./theta;