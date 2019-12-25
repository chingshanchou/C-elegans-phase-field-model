clear;

global n L h count eggphi

n=256;
L=90;
tol_count=6e3;
err=1e-9;
dx=L/n;
h=1e-3;
x0=linspace(-L/2,L/2,n+1);
x=x0(1:n);
[xx,yy]=meshgrid(x,x);
diff=1;

[phi,Bphi,my,a1,a10,a11,p,eggphi]=PF_phi_new;

k=(2*pi/L)*[0:(n/2-1) (-n/2):(-1)];
cplx=complex(0,1);                                                        % time step

[kx,ky]=meshgrid(k,k);
tt=5e-3;
Bphi=18*phi.^2.*(1-phi).^2;

Bphi0=Bphi+tt;

Bphi_k=fft2(Bphi);


Bphix=real(ifft2(cplx*Bphi_k.*kx))./Bphi0;
Bphiy=real(ifft2(cplx*Bphi_k.*ky))./Bphi0;

count=1;
kpa=0.1;
alpha=50;

colormap(jet)
surf(xx,yy,phi);
shading interp
set(gca,'FontSize',18,'FontWeight','bold');
colorbar('FontSize',18,'FontWeight','bold');
title('Cdc42-GTP');
view(2)
axis([-L/2 L/2 -L/2 L/2])
u = zeros(size(xx));
v = zeros(size(xx));


writerObj = VideoWriter('EGG.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj);

M(1)=getframe;
writeVideo(writerObj, M(1));

fileID1 = fopen('flphi.bin','a+');
fileID2 = fopen('flBphi.bin','a+');
fileID3 = fopen('fla1.bin','a+');
fileID4 = fopen('fla10.bin','a+');
fileID5 = fopen('fla11','a+');
fileID6 = fopen('flp.bin','a+');
fileID7 = fopen('flmy.bin','a+');

while (diff>err && count<=tol_count) && ~isnan(diff)
   
    [propel,phi,Bphi,diff,u,v,my,a1,a10,a11,p]=NewCelegans_my(phi,alpha,u,v,my,a1,a10,a11,p);
   
    modsize= 100;
    if mod(count,modsize)==0
        j=count/modsize+1;
        colormap(jet)
        % find the indices around the boundary 
        lambda = 10^(-3);
        ind = find(Bphi>lambda);
        my_rest = zeros(size(my));
        my_rest(ind) = my(ind);
        surf(xx,yy,my_rest);
        shading interp
        set(gca,'FontSize',18,'FontWeight','bold');
        colorbar('FontSize',18,'FontWeight','bold');
        title('Pedestal');
        view(2)
        axis([-L/2 L/2 -L/2 L/2])
        %       caxis([0,0.18])
        %M(j)=getframe;
        myFrame=getframe;
        size(myFrame.cdata)
        writeVideo(writerObj, myFrame);
        
        fwrite(fileID1,phi(:),'double');
        fwrite(fileID2,Bphi(:),'double');
        fwrite(fileID3,a1(:),'double');
        fwrite(fileID4,a10(:),'double');
        fwrite(fileID5,a11(:),'double');
        fwrite(fileID6,p(:),'double');
        fwrite(fileID7,my(:),'double');
        
    end
    count=count+1
    diff
end
%movie2avi(M,'Pedestal.avi','fps',20);
close(writerObj);

fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
fclose(fileID4);
fclose(fileID5);
fclose(fileID6);
fclose(fileID7);

%         colormap(jet)
%         surf(xx,yy,Bphi.*Am);
%         shading interp
%         set(gca,'FontSize',18,'FontWeight','bold');
%         colorbar('FontSize',18,'FontWeight','bold');
%         title('Pedestal');
%         view(2)
