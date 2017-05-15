
% code to compare approximated PolSpice like Cl determination with full
% case

%%

% load wigner small d arrays, variables:
% 'd_theta','d_ell','d_p2p2','d_m2p2'
load('wigner.mat')

%%
figure(1),clf
figure(2),clf

%read in xi's (from a default cosmology, full spherical, non-Hankel case)
m=dlmread('xipm_nonlimber_spherical_wigner.dat'); 
theta=m(:,1);
xip_full=m(:,2);
xim_full=m(:,3);
theta=theta.*300./max(theta); %apply transform for file format 

figure(1)
%make some plots 
subplot(1,2,1),...
y1=theta.*xip_full;
semilogx(theta,y1./1e-4,'k'),hold on
set(gca,'FontSize',13.),xlabel('$\theta$/arcmin','Interpreter','Latex'),ylabel('$\theta\xi_+(\theta)/10^{-4}$','Interpreter','Latex')

subplot(1,2,2),...
y1=theta.*xim_full;
semilogx(theta,y1./1e-4,'k'),hold on
set(gca,'FontSize',13.),xlabel('$\theta$/arcmin','Interpreter','Latex'),ylabel('$\theta\xi_-(\theta)/10^{-4}$','Interpreter','Latex')

ntheta=size(theta,1);
the=theta.*pi./180./60.; %convert to radians 

kp=1.;

lmin=100; 
lmax=2000; 
nl=200;
clE=zeros(nl,1);
clEa=zeros(nl,1);
lpp=zeros(nl,1);
for l=1:nl
    ell=10.^(log10(lmin)+(l-1.).*(log10(lmax)-log10(lmin))./(nl-1.));
    lpp(l)=ell;
    fprintf('%f\n',ell);
    for t=1:ntheta 
        if (t>1), dtheta=the(t)-the(t-1); else dtheta=the(t); end
        
        % interpolate tabulared wigner small-d functions to desired values
        lindex=max(find(d_ell<ell));
        thetaindex=max(find(d_theta<the(t)));
        if isempty(thetaindex)==1 
            d22=0.;
            dm22=0.;
        else
            d22=d_p2p2(lindex,thetaindex);
            dm22=d_m2p2(lindex,thetaindex);
        end
                    
        clE(l)=clE(l)+2.*pi.*(xip_full(t).*d22.*kp+xim_full(t).*dm22.*(1.-kp)).*sin(the(t)).*dtheta;  
        clEa(l)=clEa(l)+2.*pi.*(xip_full(t).*besselj(0,ell.*the(t)).*kp+xim_full(t).*besselj(4,ell.*the(t)).*(1.-kp)).*the(t).*dtheta; 
        
    end
end
    
figure(2)
subplot(1,2,1),...
loglog(lpp,lpp.*lpp.*clE.*10^5,'k'),hold on
loglog(lpp,lpp.*lpp.*clEa.*10^5,'b'),hold on
set(gca,'FontSize',13.),xlabel('$\ell$ mode','Interpreter','Latex'),ylabel('$\ell^2 C^E(\ell) 10^5$','Interpreter','Latex')
legend('Full','Approximated','Location','NorthWest') 
xlim([100 2000])
ylim([0.1 1000.])

subplot(1,2,2),...
semilogx(lpp,clE./clEa-1.,'k'),hold on
ylim([-2.5 2.5])
set(gca,'FontSize',13.),xlabel('$\ell$ mode','Interpreter','Latex'),ylabel('$C^{E, Full}(\ell)/C^{E, Approx}(\ell)-1$','Interpreter','Latex')
semilogx(lpp,zeros(size(lpp,1),1),'k')
xlim([100 2000])



