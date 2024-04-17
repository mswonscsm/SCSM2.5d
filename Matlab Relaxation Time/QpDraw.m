function Q_f=QpDraw(QPval, xfrq, tau_ep,tau_sig)
N=length(tau_ep);
K=4*N;
wmin=xfrq(1); wmax=xfrq(end);
wk=wmin*(wmax/wmin).^([0:K-1]/(K-1));
QPval_K=interp1(xfrq, QPval, wk);
theta = 1./tau_sig;
kl=(tau_ep./tau_sig-1)./N;
Q_freq1=zeros(1,size(xfrq,2));
Q_freq2=zeros(1,size(xfrq,2));

for l=1:N
    Q_freq1=Q_freq1+xfrq.*theta(l)*kl(l)./(theta(l)^2+xfrq.^2);
    Q_freq2=Q_freq2+xfrq.^2*kl(l)./(theta(l)^2+xfrq.^2);
end
QB=zeros(K,1);
QB(:)=1./QPval_K;
Q_f=Q_freq1./(1+Q_freq2);
