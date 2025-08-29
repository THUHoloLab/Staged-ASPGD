function [U_in] = Fx_Guassian(R_0)
global lambda N L0 dx0 f pz

x0=(-N/2+1:1:N/2)*dx0;
y0=(-N/2+1:1:N/2)*dx0;
[xx0,yy0]=meshgrid(x0,y0);
ap=double((xx0.^2+yy0.^2)<=(L0/2)^2);
[phi,r] = cart2pol(xx0,yy0);

% U_inx=exp(-xx0.^2/(2*R_0^2)); % R_0 is the 1/e radius
% U_iny=exp(-yy0.^2/(2*R_0^2));

U_inx=exp(-xx0.^2/(R_0^2));% R_0 is the 1/(e2) radius of the intensity
U_iny=exp(-yy0.^2/(R_0^2));
U_in=U_inx.*U_iny.*ap;
end

