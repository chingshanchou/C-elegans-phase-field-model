function [propel,pnew,Bphi,diff,u,v,my,a1,a10,a11,p]=NewCelegans_my(phi,alpha,u,v,my,a1,a10,a11,p)

global n L eps count eggphi h

% basic settings
dx=L/n;
k=(2*pi/L)*[0:(n/2-1) (-n/2):(-1)];
cplx=complex(0,1);                                                        % time step
%h=1e-3;

x=linspace(-L/2,L/2,n);
[X,Y]=meshgrid(x,x);

[kx,ky]=meshgrid(k,k);
ksq=kx.^2+ky.^2;

G1=36*(phi.*(1-phi).^2-phi.^2.*(1-phi));
G2=36*((1-phi).^2-4*phi.*(1-phi)+phi.^2);

% parameters for numerical experiments
Dm=0.01; v0=0.1; tau=2; gamma=0.4; Ma=10; A0=49.7316; theta=pi/4;
v0x = v0;
v0y =v0;

phi_k=fft2(phi);
phi_x=real(ifft2(cplx*phi_k.*kx));
phi_y=real(ifft2(cplx*phi_k.*ky));
abs_g=sqrt(phi_x.^2+phi_y.^2); % gradient phi

tt=5e-3;
Bphi=18*phi.^2.*(1-phi).^2;
Bphi_k=fft2(Bphi);
Bphi0=Bphi+tt;
phi0=phi+tt;

phix=real(ifft2(cplx*phi_k.*kx))./phi0;
phiy=real(ifft2(cplx*phi_k.*ky))./phi0;

Bphix=real(ifft2(cplx*Bphi_k.*kx))./Bphi0;
Bphiy=real(ifft2(cplx*Bphi_k.*ky))./Bphi0;

rphi=abs_g./phi0;
lap_phi = real(ifft2(-ksq.*phi_k));

[a1,a10,a11,p,my,my_x,my_y] = f_react_diff(Bphix,Bphiy,a1,a10,a11,p,my);

my_old = my;
my_oldip1 =[my_old(:,2:end),my_old(:,end-1)];
my_oldim1 =[my_old(:,2),my_old(:,1:end-1)];
my_oldjp1 =[my_old(2:end,:);my_old(end-1,:)];
my_oldjm1 =[my_old(2,:);my_old(1:end-1,:)];

h1=L/n;
my_oldx =(my_oldip1-my_oldim1)/(2*h1);
my_oldy =(my_oldjp1-my_oldjm1)/(2*h1);
%my = hill(my);

my_grad = sqrt(my_x.^2+my_y.^2);

myip1 =[my(:,2:end),my(:,end-1)];
myim1 =[my(:,2),my(:,1:end-1)];
myjp1 =[my(2:end,:);my(end-1,:)];
myjm1 =[my(2,:);my(1:end-1,:)];

h1=L/n;
my_x =(myip1-myim1)/(2*h1);
my_y =(myjp1-myjm1)/(2*h1);

[u,v]= u_velocity(phi,Bphi,phi_x,phi_y,my,my_old,my_x,my_y,my_oldx,my_oldy,kx,ky,u,v,G1,lap_phi,p);
my = my_old;
alpha = 1;
% cell area

%%%%%%%%%%%%%%%%
A=sum(sum(abs(phi).*dx^2));
Ln=sum(sum(((phi_x.^2+phi_y.^2) +G1/eps^2).*dx^2));
A0 = 523.7885;
L0 = 51.3693;
Ma = 20;
Ml =10;
area=fft2(-Ma*(A-A0).*abs_g);
length =fft2(-Ml*(Ln-L0).*abs_g);
%%%%%%%%%%%%%%%%

tension=gamma*fft2(-G1/eps);
propel=fft2(alpha.*(-phi_x.*u-phi_y.*v));

cur_term = gamma*fft2(sqrt(eps)*(lap_phi.*eps.*abs_g-G1/eps.*abs_g));
eggcoeff = 20;
%eggcoeff = 0;

egg_term = -eggcoeff*fft2(phi.*eggphi);
ptmp=(phi_k+h*(tension+propel+cur_term))./(1+h*eps*gamma*ksq);

pnew=real(ifft2(ptmp));

% error calculating
diffp=sqrt(dx^2*sum(sum(abs(pnew-phi).^2)));
diff=diffp;

end
