function [phi1,Bphi,my,a1,a10,a11,p,eggphi]=PF_phi_new

global n L eps r

r=1;
eps=L/n;
eps=2;
x=linspace(-L/2,L/2,n);
[X,Y]=meshgrid(x,x);

eggdiff= 5.5;

rx = 15;
ry = 25;

smootstep = @(x) tanh(x/(0.1*(eps/4)))+1;
phi1 = 0.5*smootstep(1-sqrt((X/(rx)).^2+((Y)/(ry)).^2));



eggphi = -0.5*smootstep(1-sqrt((X/(rx+eggdiff)).^2+((Y)/(ry+eggdiff)).^2))+1;

a1high = 0.5162;
a1low = 0.0164;
a10high = 0.2812;
a10low = 0.0086;
a11high = 0.4971;
a11low = 0.0069;
plow = 0.0525;
phigh = 0.6532;
mycons = 0.7;

myhigh = 0.7;
mylow = 0.2;
Kin = ry-5;

div = 2;
a1 = a1high*(1-(tanh((-(Y+Kin))/(eps/div))+1)/2)+a1low;
a10 = a10high*(1-(tanh((-(Y+Kin))/(eps/div))+1)/2)+a10low;
a11 = a11high*(1-(tanh((-(Y+Kin))/(eps/div))+1)/2)+a11low;
my = myhigh*(1-(tanh((-(Y+Kin))/(eps/div))+1)/2)+mylow;

p = phigh*(tanh((-(Y+Kin))/(eps/div))+1)/2+plow;

Bphi=18*phi1.^2.*(1-phi1).^2;
end
