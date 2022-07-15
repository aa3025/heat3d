% Heat equation in 3D, explicit central difference in space, 
% steady state and transient in time with forward Euler method.

clear; close all; clc;


%% material props
material='copper';

switch material
    case 'water' %  very bad conductor of heat
    
    k = 0.6;                                         % thermal conductivity (W/m.K)
    cp = 4.187e+3;                                   % specific heat (J/kg K)
    rho = 1000;                                      % density (kg/m^3)

    case 'copper' % very good conductor of heat
    k = 401;                                        % thermal conductivity (W/m.K)
	cp = 390;                                       % specific heat (J/kg K)
	rho = 8940;                                     % density (kg/m^3) 
end
alph = k/(cp*rho);                                 % alpha (1/sec) (heat diffusivity)


video=0; % record video? 
videofilename='heat_eq_3D.mp4';
videoformat='MPEG-4';
if video>0
    writerObj = VideoWriter(videofilename, videoformat);  % (see also lines 179+)
    open(writerObj); 
end

%% domain
Lx= 1;                  % plate width (m)
Ly= 1;                  % plate depth (m)
Lz= 1;                  % plate height (m)
Nx=51;                 % nodes in x direction (odd number to have central plane to slice through)
Ny=51;                 % nodes in y direction (odd number to have central plane to slice through)
Nz=51;                  % nodes in z direction (odd number to have central plane to slice through)
dx=Lx/Nx;               % spacing along x
dy=Ly/Ny;               % spacing along y
dz=Lz/Nz;               % spacing along z

xc=(Nx+1)/2;% center plane of the domain
yc=(Ny+1)/2;
zc=(Nz+1)/2;

dt = 1/(2*alph*((1/dx^2)+(1/dy^2)+(1/dz^2)));	% stable time step, sec
tmax=3600; % max time, sec

%% Initial and boundary conditions (1st type, a.k.a Dirichlet)
T_0= 0;                 % Initial temperature in all nodes 
T_right   = -20;        % temperature at x=Lx
T_left   = 20;        % temperature at x=0 
T_top  = 0;             % temperature at y=0 
T_bott  = 0;            % temperature at y=Ly  

tol_ss=0.01;            % tolerance for steady state

%% Constants

iter=1;  % counter of iterations
error  = 1; % initial error for iterations

%% Initial conditions for Steady State
Tss=zeros(Ny,Nx,Nz);

Tss(:,1,:)    =T_left;
Tss(:,Nx,:)   =T_right;

Tss(1,:,:)    =T_bott;
Tss(Ny,:,:)   =T_top;

%% plotting 
x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
z=linspace(0,Lz,Nz);
[X Y Z]=meshgrid(x,y,z);

figure('units','pixels','position',[100 100 1280 720])

 
%% Steady-State Solution
tic
i=2:Ny-1; %inner nodes along y
j=2:Nx-1; %inner nodes along x
k=2:Nx-1; %inner nodes along z
Tss_new=Tss;
T_min=T_right;
T_max=T_left;
numisosurf=25; % number of isosurfaces
isovalues=linspace(T_min,T_max,numisosurf);
tmax=3600; % max time, sec 
t=0;
k1=alph*dt; 
method='steadystate';

while t<=tmax%error >= tol_ss % by tmax or by tolerance (for Staeady State)
   
  switch method
        case 'steadystate'
    
   % Solving Steady-State Laplace's eqn on 3D FD stencil
    Tss_new(i,j,k)=(dz^2*(Tss(i,j,k-1)+Tss(i,j,k+1))+dy^2*(Tss(i+1,j,k)+Tss(i-1,j,k))+dx^2*(Tss(i,j+1,k)+Tss(i,j-1,k)))/(2*(dx^2+dy^2+dz^2));    
  
        case 'Euler'
    
    % FTCS transient -- forward Euler
   Tss_new(i,j,k) =Tss(i,j,k)+k1*(((Tss(i-1,j,k)-2*Tss(i,j,k)+Tss(i+1,j,k))/dy^2)+((Tss(i,j-1,k)-2*Tss(i,j,k)+Tss(i,j+1,k))/dx^2)+((Tss(i,j,k-1)-2*Tss(i,j,k)+Tss(i,j,k+1))/dz^2));
   
  end
    iter=iter+1;
    t=iter*dt; % pointless for Steady State
    error=max(max(max(abs(Tss_new-Tss))));    
    Tss=Tss_new;

if mod(iter,5)==0
    plot3D(X,Y,Z,Lx,Ly,Lz,Tss,isovalues,iter,t); % e.g. plot every 5 iter
end

if video >0 
         frame = getframe(gcf);
         writeVideo(writerObj,frame);
     end

 end
 
%% 
toc 
% plot3D(X,Y,Z,Lx,Ly,Lz,Tss,isovalues,iter); % final plot

 if video >0  close(writerObj); end
 


disp('COMPLETE!');
























