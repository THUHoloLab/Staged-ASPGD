function Zern_N= Fx_Zernike(R_0)
%FX_ZERNIKE generate Zernike Polynomials
% n=1 m= -1 1 (s=1,2)
% n=2 m=-2 0 2 (s=3,4,5)
% n=3 m=-3 -1 1 3 (s=6,7,8,9)
% n=4 m=-4 -2 0 2 4 (s=10 11 12 13 14)
% n=5 m=-5 -3 -1 1 3 5 (s=15 16 17 18 19 20)
% n=6 m=-6 -4 -2 0 2 4 6 (s=21 22 23 24 25 26 27 )
% n=7 m=-7 -5 -3 -1 1 3 5 7 (s=28 29 30 31 32 33 34 35)

% R_0 is the radius of input Guassian beam

global N dx0
M=ceil(2*R_0/dx0);
if mod((N-M),2)~=0
    M=M+1;
end
x=linspace(-1,1,M); % normalized radius, from -1 to 1 
y=linspace(1,-1,M);
[x,y]=meshgrid(x,y);
[t,r] = cart2pol(x,y);
zern_mat=zeros([M,M,35]);
ap=double(r<=1);

s=1;
for n = 1:15
    for m = -n:2:n            
            Z= zernike(r,t,n,m).*ap;
            zern_mat(:,:,s)=Z/max(abs(Z),[],'all');
                % figure, imagesc(elliptical_crop(zern_mat(:,:,s),1));
                % title(strcat('n = ',{' '},string(n),', m = ',{' '},string(m),'s = ',{' '},string(s)));
                % colormap jet;
                % drawnow;
     s=s+1;   
     end
end

Zern_N=padarray(zern_mat,[(N-M)/2 (N-M)/2]);
% Zern_N=zern_mat;

end

