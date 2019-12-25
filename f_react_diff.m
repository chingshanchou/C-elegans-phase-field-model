function [a1,a10,a11,p,my,my_x,my_y] = f_react_diff(Bphix,Bphiy,a1,a10,a11,p,my)

global L n count

%dt=1e-2;
dt=2e-1;
h=L/n;

alpha_y = 1;
pho_y = 1;

beta_1 = 0.011979;
beta_2 = 4.6526;
beta_3 = 4.177344;
beta_4 = 0.911258;
beta_5 = 0.011441;
beta_6 = 5.003463;
beta_7 = 0.073918;
beta_8 = 0.082689;
beta_9 = 0.724544;
beta_10 = 0.1;
beta_11 = 0.1;
beta_12 = 0.1;

mu = 0.002;
D1 = 0.002;
D2 = 0.0015;
D3 = 0.002;

myip1 =[my(:,2:end),my(:,end-1)];
myim1 =[my(:,2),my(:,1:end-1)];
myjp1 =[my(2:end,:);my(end-1,:)];
myjm1 =[my(2,:);my(1:end-1,:)];

my_x =(myip1-myim1)/(2*h);
my_y =(myjp1-myjm1)/(2*h);

a1ip1 =[a1(:,2:end),a1(:,end-1)];
a1im1 =[a1(:,2),a1(:,1:end-1)];
a1jp1 =[a1(2:end,:);a1(end-1,:)];
a1jm1 =[a1(2,:);a1(1:end-1,:)];

a1_x =(a1ip1-a1im1)/(2*h);
a1_y =(a1jp1-a1jm1)/(2*h);

a10ip1 =[a10(:,2:end),a10(:,end-1)];
a10im1 =[a10(:,2),a10(:,1:end-1)];
a10jp1 =[a10(2:end,:);a10(end-1,:)];
a10jm1 =[a10(2,:);a10(1:end-1,:)];

a10_x =(a10ip1-a10im1)/(2*h);
a10_y =(a10jp1-a10jm1)/(2*h);

a11ip1 =[a11(:,2:end),a11(:,end-1)];
a11im1 =[a11(:,2),a11(:,1:end-1)];
a11jp1 =[a11(2:end,:);a11(end-1,:)];
a11jm1 =[a11(2,:);a11(1:end-1,:)];

a11_x =(a11ip1-a11im1)/(2*h);
a11_y =(a11jp1-a11jm1)/(2*h);

pip1 =[p(:,2:end),p(:,end-1)];
pim1 =[p(:,2),p(:,1:end-1)];
pjp1 =[p(2:end,:);p(end-1,:)];
pjm1 =[p(2,:);p(1:end-1,:)];

p_x =(pip1-pim1)/(2*h);
p_y =(pjp1-pjm1)/(2*h);

h_inv =1/h;
a1xfor = (a1ip1-a1)*h_inv;
a1xback = (a1-a1im1)*h_inv;

a1yfor = (a1jp1-a1)*h_inv;
a1yback = (a1-a1jm1)*h_inv;

a10xfor = (a10ip1-a10)*h_inv;
a10xback = (a10-a10im1)*h_inv;

a10yfor = (a10jp1-a10)*h_inv;
a10yback = (a10-a10jm1)*h_inv;

a11xfor = (a11ip1-a11)*h_inv;
a11xback = (a11-a11im1)*h_inv;

a11yfor = (a11jp1-a11)*h_inv;
a11yback = (a11-a11jm1)*h_inv;

pxfor = (pip1-p)*h_inv;
pxback = (p-pim1)*h_inv;

pyfor = (pjp1-p)*h_inv;
pyback = (p-pjm1)*h_inv;

mxfor = (myip1-my)*h_inv;
mxback = (my-myim1)*h_inv;
mxcen = (myip1-myim1)/(h^2);

myfor = (myjp1-my)*h_inv;
myback = (my-myjm1)*h_inv;
mycen = (myjp1-myjm1)/(h^2);

mxpos = (mxcen>0).*ones(size(my));
mxneg = (mxcen<0).*ones(size(my));
mypos = (mycen>0).*ones(size(my));
myneg = (mycen<0).*ones(size(my));

mxposcen = mxpos.*mxcen;
mxnegcen = mxneg.*mxcen;
myposcen = mypos.*mycen;
mynegcen = myneg.*mycen;

grada1gradm =  mxposcen.*a1xback+mxnegcen.*a1xfor+myposcen.*a1yback+mynegcen.*a1yfor;
grada10gradm =  mxposcen.*a10xback+mxnegcen.*a10xfor+myposcen.*a10yback+mynegcen.*a10yfor;
grada11gradm =  mxposcen.*a11xback+mxnegcen.*a11xfor+myposcen.*a11yback+mynegcen.*a11yfor;
gradpgradm = mxposcen.*pxback+mxnegcen.*pxfor+myposcen.*pyback+mynegcen.*pyfor;
gradmgradm = mxposcen.*mxback+mxnegcen.*mxfor+myposcen.*myback+mynegcen.*myfor;
% Laplacian

Lap_m = (myip1+myim1+myjp1+myjm1-4*my)/(h^2);
Lap_a1 = (a1ip1+a1im1+a1jp1+a1jm1-4*a1)/(h^2);
Lap_a10 = (a10ip1+a10im1+a10jp1+a10jm1-4*a10)/(h^2);
Lap_a11 = (a11ip1+a11im1+a11jp1+a11jm1-4*a11)/(h^2);
Lap_p = (pip1+pim1+pjp1+pjm1-4*p)/(h^2);

%% diffusion

a1_diff = D1*(Bphix.*a1_x+Bphiy.*a1_y + Lap_a1);
a10_diff = D1*(Bphix.*a10_x+Bphiy.*a10_y + Lap_a10);
a11_diff = D1*(Bphix.*a11_x+Bphiy.*a11_y + Lap_a11);
my_diff = D2*(Bphix.*my_x+Bphiy.*my_y + Lap_m);
p_diff = D3*(Bphix.*p_x+Bphiy.*p_y + Lap_p);
%% advection:

a1_adv = -mu*(a1.*(Bphix.*my_x+Bphiy.*my_y) + grada1gradm+a1.*Lap_m);
a10_adv = -mu*(a10.*(Bphix.*my_x+Bphiy.*my_y) + grada10gradm+a10.*Lap_m);
a11_adv = -mu*(a11.*(Bphix.*my_x+Bphiy.*my_y) + grada11gradm+a11.*Lap_m);
p_adv = -mu*(p.*(Bphix.*my_x+Bphiy.*my_y) + gradpgradm+p.*Lap_m);
my_adv = -mu*(my.*(Bphix.*my_x+Bphiy.*my_y) + gradmgradm+my.*Lap_m);
%%

a1 = a1+ dt*(beta_1*alpha_y - a1 -2*beta_2*a1.^2 + 2*beta_3*a11 - beta_2*alpha_y*a1 + beta_3*a10 - beta_4*p.*a1 + a1_diff + a1_adv);

a10 = a10+ dt*(beta_5*alpha_y^2 - a10 + beta_2*alpha_y*a1 -beta_3*a10 + a11- beta_6*a10  - beta_4*p.*a10+ a10_diff + a10_adv);

a11 = a11+ dt*(beta_2*a1.^2 -beta_3*a11 -a11 + beta_6*a10-2*beta_4*p.*a11+ a11_diff + a11_adv);

p = p+ dt*(beta_7*pho_y -beta_8*p- beta_9*(a1 +a10 +2*a11 ).*p + p_diff + p_adv);

my = my+ dt*(beta_10.*(beta_11./(beta_11+p))  -beta_12*my+ + my_diff + my_adv);


end