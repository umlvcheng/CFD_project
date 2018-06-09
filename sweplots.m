function sweplots(u,v,eta,xc,yc,xe,ye,n,dt,dx,dy)

figure(1)
set(gcf,'Position',[50 50 800 1200])

subplot 311
%surf(1:M,1:N,ones([M N]),eta');zlim([-0.1 0.1])
contourf(xc,yc,eta');%zlim([-0.1 0.1])
grid off;shading interp;
c = (1:256)'/256;
cyan = [0*c c.^2 c.*(2-c)];
colormap(cyan);colorbar;caxis([-0.2 0.2]);%title(colorbar,'(m)')
%title(['\eta, t =' num2str(n*dt) 's, dx = ' num2str(dx) ', dy = ' num2str(dy)]);
%xlabel 'length [ m]'; ylabel 'width [ m]';
set(gca,'Fontsize',20,'FontName','Cambria');
    

subplot 312
contourf(xe,yc,u',10);%zlim([0.05 -0.05]);
grid off;shading interp;
c = (1:512)'/512;
cyan = [0*c c.^2 c.*(2-c)];
colormap jet;colorbar;caxis([-0.35 0.2]);%title(colorbar,'(m/s)')
title 'u'
%colorbar;colormap jet;%caxis([-0.05 0.05])
%xlabel 'length [ m]'; ylabel 'width [ m]';
set(gca,'Fontsize',20,'FontName','Cambria');

subplot 313
contourf(xc,ye,v',10);grid off;shading interp;
c = (1:512)'/512;
cyan = [0*c c.^2 c.*(2-c)];
colormap(jet);colorbar;caxis([-0.08 0.08]);%title(colorbar,'(m/s)')
title 'v'
%xlabel 'length [ m]'; ylabel 'width [ m]';
set(gca,'Fontsize',20,'FontName','Cambria');


drawnow
