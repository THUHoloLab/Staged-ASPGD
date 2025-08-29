clear;
clc;
close all;
%% Model settings
%mm
global lambda N L0 dx0 f
lambda=1064e-6;

NA=0.45;
f=10;

R=NA*f;
L0=2*R;

dx0=12.5e-3;
dx1=lambda*f/L0;

N=ceil(L0/dx0);

L_flat=0.1;% half length
R_0=3/2;

% Zernike modes
Zern_N=Fx_Zernike(R);

% axis
x0=(-N/2+1:1:N/2)*dx0;
y0=(-N/2+1:1:N/2)*dx0;
[xx0,yy0]=meshgrid(x0,y0);
ap=double((xx0.^2+yy0.^2)<=(L0/2)^2);
[phi,r] = cart2pol(xx0,yy0);

x1=(-N/2+1:1:N/2)*dx1;
y1=(-N/2+1:1:N/2)*dx1;
[xx1,yy1]=meshgrid(x1,y1);
% Uin
U_in= Fx_Guassian(R_0);
figure;imshow(ap.*U_in.^2,[]);title('Iin');
%% phase abbration
Zern_AP=12;
Noise_in=pi*[0 0 0.2 -0.2 0.3 0.2 0.15 -0.2 0.1 0.15 0.05 -0.15];Noise_in=reshape(Noise_in,1,1,Zern_AP);
phase_abb=sum(pi*Noise_in.*Zern_N(:,:,1:Zern_AP),3);
%% amplitude abbration
U_inx=exp(-xx0.^2/(1.2*R_0^2));% R_0 is the 1/(e2) radius of the intensity
U_iny=exp(-yy0.^2/(0.8*R_0^2));
E_abb=U_inx.*U_iny.*ap;
figure;
subplot(1,2,1);imagesc(E_abb);title('E\_abb');axis image;axis off;
subplot(1,2,2);imagesc(U_in);title('U\_in');axis image;axis off;
%% Target
I_t=exp(-((xx1/L_flat).^(20)+(yy1/L_flat).^(20)));
%% Phase for Ideal flattop beam
% this part used the 1/(e2) radius
beta=2*sqrt(2*pi)*R_0*L_flat/(f*lambda);
ens_x=sqrt(2)*x0/R_0;
ens_y=sqrt(2)*y0/R_0;

% phi_x0 1XN
phi_x0=1/erf(max(abs(ens_x)))*(beta)*...
   ((sqrt(pi)/2).*ens_x.*erf(ens_x)...
   +1/2*exp(-(ens_x).^2)...
   -1/2);

% phi_y0 1XN
phi_y0=1/erf(max(abs(ens_y)))*(beta)*...
   ((sqrt(pi)/2).*ens_y.*erf(ens_y)...
   +1/2*exp(-(ens_y).^2)...
   -1/2);

[phi_xx0,phi_yy0]=meshgrid(phi_x0,phi_y0);

phi_in=phi_xx0+phi_yy0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

I_ideal=abs(fftshift(fft2(ifftshift(U_in.*exp(1j*phi_in))))).^2;
figure;
subplot(2,2,1);imshow(mod(phi_in,2*pi),[]);title('Stationary phase');
subplot(2,2,2);imshow(I_ideal,[]);title('Output');
subplot(2,2,3);imshow(U_in,[]);title('Input');
subplot(2,2,4);plot(x1,I_ideal(N/2+1,:));xlabel('x1/mm');

%% Draw pictures before optimization
E_in = ap.*E_abb.*exp(1j*phase_abb); 

% PSF with abbration
PSF_A=abs(fftshift(fft2(ifftshift(E_in)))).^2;
figure;imshow(PSF_A,[]);title('PSF_A');

SLM_Pha=phi_in;
I_o=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*SLM_Pha))))).^2;

figure;
subplot(2,2,1);
imagesc(x0,y0,(abs(E_in)).^2);title(['I\_in r(1/e2)=',num2str(R_0),'mm']);colormap("gray");xlabel('x/mm');ylabel('y/mm');axis image;
title('Input Intensity');
subplot(2,2,2);
imagesc(x0,y0,mod(phase_abb,2*pi)/(2*pi));title(['I\_in r(1/e2)=',num2str(R_0),'mm']);colormap("gray");xlabel('x/mm');ylabel('y/mm');axis image;
title('Input Phase');
subplot(2,2,3);
imagesc(x1,y1,I_o);title(['I\_in r(1/e2)=',num2str(R_0),'mm']);colormap("gray");xlabel('x/mm');ylabel('y/mm');axis image;
title('Output');
subplot(2,2,4);
plot(x1,I_o(N/2,:));
drawnow;
%% ASPGD  
tic
% initialize
Zern_Level=65;%135

beta_1=0.2;
beta_2=0.999;

enss=1e-20;

gama=100;

% OptimizeStep 1: PSF 2: Fx_MSE
for OptimizeStep=[1 2]% First optimize PSF

if  OptimizeStep==1%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++PSF
Iteration_N=200;
Coefficient_Zernike=zeros(1,1,Zern_Level);
P_0=zeros(N);
I_t=zeros(N);

I_t(N/2+1,N/2+1)=1;
I_t=I_t*sum(PSF_A,'all');

% CostFunc
CostV=zeros(1,Iteration_N);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*P_0))))).^2;
CostV(1)=Fx_STD(F3d);

elseif OptimizeStep==2
Iteration_N=2000;
P_0=phi_in;

Coefficient_Zernike=Cof_PSF;

Result_Pha=P_0+sum(Coefficient_Zernike.*Zern_N(:,:,1:Zern_Level),3);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*Result_Pha))))).^2;

I_t=exp(-((xx1/L_flat).^(20)+(yy1/L_flat).^(20)));
I_t=I_t*sum(PSF_A,'all')/sum(I_t,'all');
% CostFunc
CostVM=zeros(1,Iteration_N);
Cost_t=Fx_MSE(F3d,I_t);
alpha=1/Cost_t;
CostVM(1)=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);
end

%stepset
if OptimizeStep==1
Noiseamplitude=0.1;
% estimate the Steplength
RandomNoise=Noiseamplitude*(-1).^binornd(1,0.5,1,1,Zern_Level);
AddCoeff=Coefficient_Zernike+RandomNoise;
AddPha=P_0+sum(AddCoeff.*Zern_N(:,:,1:Zern_Level),3);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*AddPha))))).^2;

AddCostV=Fx_STD(F3d);

SubCoeff=Coefficient_Zernike-RandomNoise;
SubPha=P_0+sum(SubCoeff.*Zern_N(:,:,1:Zern_Level),3);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*SubPha))))).^2;

SubCostV=Fx_STD(F3d);

delta_J=AddCostV-SubCostV;
g=delta_J*RandomNoise/(2*Noiseamplitude^2);%g 1*1*level

veloc=zeros(1,1,Zern_Level);
momentum=zeros(1,1,Zern_Level);

momentum=beta_1*momentum+(1-beta_1)*g;
veloc=beta_2*(veloc)+(1-beta_2)*(g.^2);
momentum_v=momentum/(1-beta_1^(1));
veloc_v=veloc/(1-beta_2^(1));

deta_CostV=momentum./sqrt(veloc+enss);
StepLength=10^(-2);% + maxCostV, - min
                       
elseif OptimizeStep==2
Noiseamplitude=0.3;%1
% estimate the Steplength
RandomNoise=Noiseamplitude*(-1).^binornd(1,0.5,1,1,Zern_Level);
AddCoeff=Coefficient_Zernike+RandomNoise;
AddPha=P_0+sum(AddCoeff.*Zern_N(:,:,1:Zern_Level),3);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*AddPha))))).^2;

AddCostV=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);

SubCoeff=Coefficient_Zernike-RandomNoise;
SubPha=P_0+sum(SubCoeff.*Zern_N(:,:,1:Zern_Level),3);
F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*SubPha))))).^2;

SubCostV=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);

delta_J=AddCostV-SubCostV;
g=delta_J*RandomNoise/(2*Noiseamplitude^2);%g 1*1*level

veloc=zeros(1,1,Zern_Level);
momentum=zeros(1,1,Zern_Level);

momentum=beta_1*momentum+(1-beta_1)*g;
veloc=beta_2*(veloc)+(1-beta_2)*(g.^2);
momentum_v=momentum/(1-beta_1^(1));
veloc_v=veloc/(1-beta_2^(1));

% exponential
deta_CostV=momentum_v./sqrt(veloc_v+enss);
Min_gradient=min(abs(deta_CostV),[],'all');
StepLength=-10^(-2);
end

veloc=zeros(1,1,Zern_Level);
momentum=zeros(1,1,Zern_Level);
for ii=1:Iteration_N
    RandomNoise=Noiseamplitude*(-1).^binornd(1,0.5,1,1,Zern_Level);

    AddCoeff=Coefficient_Zernike+RandomNoise;
    AddPha=P_0+sum(AddCoeff.*Zern_N(:,:,1:Zern_Level),3);
    F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*AddPha))))).^2;
    
    if OptimizeStep==1   
    AddCostV=Fx_STD(F3d);
    elseif OptimizeStep==2   

    AddCostV=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);
    end
    
    SubCoeff=Coefficient_Zernike-RandomNoise;
    SubPha=P_0+sum(SubCoeff.*Zern_N(:,:,1:Zern_Level),3);
    F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*SubPha))))).^2;
    
    if OptimizeStep==1   
    SubCostV=Fx_STD(F3d);
    elseif OptimizeStep==2
 
    SubCostV=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);
    end

    delta_J=AddCostV-SubCostV;
    g=delta_J*RandomNoise/(2*Noiseamplitude^2);

    momentum=beta_1*momentum+(1-beta_1)*g;
    veloc=beta_2*(veloc)+(1-beta_2)*(g.^2);
    momentum_v=momentum/(1-beta_1^(ii));
    veloc_v=veloc/(1-beta_2^(ii));

    deta_CostV=momentum_v./sqrt(veloc_v+enss); 
           
    Coefficient_Zernike= Coefficient_Zernike + StepLength.*deta_CostV;        
      
    Result_Pha=P_0+sum(Coefficient_Zernike.*Zern_N(:,:,1:Zern_Level),3);
    F3d=abs(fftshift(fft2(ifftshift(E_in.*exp(1j*Result_Pha))))).^2;
    
    
    if OptimizeStep==1
    CostV(ii+1)=Fx_STD(F3d);
    CostV(ii+1)
    elseif OptimizeStep==2
    CostVM(ii+1)=alpha*Fx_MSE(F3d,I_t)+gama*(1-Fx_Structure(F3d,I_t)).^(0.1);
    CostVM(ii+1)
    end
    
end
if OptimizeStep==1
    Cof_PSF=Coefficient_Zernike;
elseif OptimizeStep==2
    Cof_MSE=Coefficient_Zernike;
end

figure;
subplot(2,3,1);imshow(I_t,[]);title('Target');
if  OptimizeStep==2
subplot(2,3,2);imagesc(x1,y1,F3d);title('After Optimization');colormap gray;axis image;axis off;
line([x1(1);x1(end)],[y1(N/2+1);y1(N/2+1)],[max(F3d,[],'all');max(F3d,[],'all')],'LineWidth',1.2,'Color','y','LineStyle','--');
subplot(2,3,3);imagesc(x1,y1,I_o);title(['Before Optimization, U=']);colormap gray;axis image;axis off; 
line([x1(1);x1(end)],[y1(N/2+1);y1(N/2+1)],[max(I_o,[],'all');max(F3d,[],'all')],'LineWidth',1.2,'Color','y','LineStyle','--');
subplot(2,3,5);plot(x1,F3d(N/2+1,:));title('After Optimization');xlabel('x/mm');
subplot(2,3,6);plot(x1,I_o(N/2+1,:));title('Before Optimization');xlabel('x/mm');
elseif OptimizeStep==1
subplot(2,3,2);imagesc(x1,y1,F3d);title('After Optimization');colormap gray;axis image;axis off;
subplot(2,3,3);imshow(PSF_A,[]);title('Before Optimization');   
subplot(2,3,6);imshow(mod(phase_abb,2*pi),[]);title('Phase aberration');
end
if OptimizeStep==1
subplot(2,3,4);plot(CostV(CostV~=0));title('Costvalue');
elseif OptimizeStep==2
subplot(2,3,4);plot(CostVM(CostVM~=0));title('Costvalue');
end
end
toc

