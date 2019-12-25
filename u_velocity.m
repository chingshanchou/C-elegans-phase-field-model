function [u,v]= u_velocity(phi,Bphi,phi_x,phi_y,my,my_old,my_x,my_y,my_oldx,my_oldy,kx,ky,u,v,G1,lap_phi,p)

global n L eps r count eggphi

nu0 = 1000;
cyt_vis = 10;
phi_c = cyt_vis/nu0;
%phi_tilde =2;
phi_tilde = 8;
psi = 0.5;
ct=600;
Ma = 100;
Ml=0;
cm = 50;
cg = 500;
kappa = 20;
gamma2 = 20;
h1=L/n;
h=L/n;
tt = 5e-3;


ksq=kx.^2+ky.^2;
cplx=complex(0,1);
mi = 10;
phi_k = fft2(phi);
G2=36*((1-phi).^2-4*phi.*(1-phi)+phi.^2);

%p
pip1 =[p(:,2:end),p(:,end-1)];
pim1 =[p(:,2),p(:,1:end-1)];
pjp1 =[p(2:end,:);p(end-1,:)];
pjm1 =[p(2,:);p(1:end-1,:)];

p_x =(pip1-pim1)/(2*h1);
p_y =(pjp1-pjm1)/(2*h1);

Bphi = 4*phi.*(1-phi);
%Bphi_k = fft2(Bphi);
Bphi_x=4*(phi_x-2*phi.*phi_x);
Bphi_y=4*(phi_y-2*phi.*phi_y);

abs_g = sqrt(phi_x.^2+phi_y.^2);
abs_gsq = (phi_x.^2+phi_y.^2+tt);
%normal vectors
n1 = (phi_x./abs_gsq).*(abs_g);
n2 = (phi_y./abs_gsq).*(abs_g);

% tension
tension = eps*lap_phi-G1/eps;
F_ten_u = -(gamma2+ct*my).*(tension).*phi_x;
F_ten_v = -(gamma2+ct*my).*(tension).*phi_y;

% bending
% bend=-kpa*(ksq.*fft2(G1/eps^2)-fft2(real(ifft2(-ksq.*phi_k)).*G2/eps^2)+fft2(G2.*G1/eps^4));
bend = kappa*eps*(real(ifft2( ksq.^2.*phi_k)) +real(ifft2( ksq.*fft2(G1/eps^2)))-lap_phi.*G2/eps^2 +G2.*G1/eps^4);
F_bend_u = bend.*phi_x;
F_bend_v = bend.*phi_y;

% eggshell
eggcoeff = 20000;
%eggcoeff = 0;
egg_term = eggcoeff*phi.*eggphi;

A0 =  1176;
L0 = 60.6232;

A=sum(sum(abs(phi).*h^2));
Ln=sum(sum(((phi_x.^2+phi_y.^2) +G1/eps^2).*h^2));
area=Ma*(A-A0);
length =Ml*(Ln-L0);

penalthy = (egg_term+area+length);

penalthy_u = penalthy.*n1;
penalthy_v = penalthy.*n2;

my_grad = sqrt(my_x.^2+my_y.^2);
% alignment force
alig_u = (cm*my+cg*my_grad).*n1;
alig_v = (cm*my+cg*my_grad).*n2;

% total membrane force
F_mem_u = F_ten_u+F_bend_u+penalthy_u+alig_u;
F_mem_v = F_ten_v+F_bend_v+penalthy_v+alig_v;

phi_old = phi;
phix_old = phi_x;
phiy_old = phi_y;

u_temp=u;
v_temp=v;
i=0;

while (max([u_temp(:);v_temp(:)]-[u(:);v(:)])>0.01*max([u_temp(:);v_temp(:)])|i<10)
    i=i+1;
    
    u= u_temp;
    v= v_temp;
    
    u_k = fft2(u);
    v_k = fft2(v);
    
    u_x = real(ifft2(cplx*u_k.*kx));
    u_y = real(ifft2(cplx*u_k.*ky));
    v_x = real(ifft2(cplx*v_k.*kx));
    v_y = real(ifft2(cplx*v_k.*ky));
    u_xy = real(ifft2(-u_k.*kx.*ky));
    v_xy = real(ifft2(-v_k.*kx.*ky));
    u_xx = real(ifft2(-u_k.*kx.^2));
    u_yy = real(ifft2(-u_k.*ky.^2));
    v_xx = real(ifft2(-v_k.*kx.^2));
    v_yy = real(ifft2(-v_k.*ky.^2));
    
    
    phi_c2 =4;
    % new phi (the domain of viscosity) is (4eta_m*phi*(1-phi)+eta_c*phi)/eta_m
    % which is: phi.*phi_c2*(1-phi).*my_old+ phi.*phi_c)=phi.*(phi_c2*(1-phi).*my_old+phi_c)
    % eta_c = cyt_vis,
    % phi_c = cyt_vis/nu0
    % eta_m = nu0
    
    phi_x = phi_x.*(phi_c2*(1-phi).*my_old+phi_c) + phi_c2*phi.*(-phi_x.*my_old+(1-phi).*my_oldx );
    phi_y = phi_y.*(phi_c2*(1-phi).*my_old+phi_c) + phi_c2*phi.*(-phi_y.*my_old+(1-phi).*my_oldy );
    phi = phi.*(phi_c2*(1-phi).*my_old+phi_c);
    
    %phi_tilde = 4*phi_tilde;
    vis_u = nu0*(u_xx.*(2*phi-phi_tilde) + 2*u_x.*phi_x+ u_yy.*(phi-phi_tilde)+ u_y.*phi_y +v_xy.*phi +v_x.*phi_y);
    vis_v = nu0*( v_xx.*(phi-phi_tilde)+ v_x.*phi_x +u_xy.*phi +u_y.*phi_x +v_yy.*(2*phi-phi_tilde) + 2*v_y.*phi_y);
    
    fft2RHS_u = fft2(vis_u + F_mem_u);
    fft2RHS_v = fft2(vis_v + F_mem_v);
    
    u_temp = real(ifft2(fft2RHS_u./(psi + nu0*phi_tilde*ksq)));
    v_temp = real(ifft2(fft2RHS_v./(psi + nu0*phi_tilde*ksq)));
    
    phi = phi_old;
    phi_x=phix_old;
    phi_y=phiy_old;
end

u= u_temp;
v= v_temp;

end