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
Lx= 2;                  %  width (m)
Ly= 1;                  %  depth (m)
Lz= 1;                  %  height (m)
Nx=41;                 % nodes in x direction (odd number to have central plane to slice through)
Ny=21;                 % nodes in y direction (odd number to have central plane to slice through)
Nz=21;                  % nodes in z direction (odd number to have central plane to slice through)
dx=Lx/Nx;               % spacing along x
dy=Ly/Ny;               % spacing along y
dz=Lz/Nz;               % spacing along z

xc=(Nx+1)/2;% center plane of the domain
yc=(Ny+1)/2;
zc=(Nz+1)/2;

dt = 1/(2*alph*((1/dx^2)+(1/dy^2)+(1/dz^2)));	% stable time step, sec for Euler ethod
tmax=3600; % max time, sec

%% Initial and boundary conditions (1st type, a.k.a Dirichlet)
% T(y,x,z) % structure of 3D array

T_0= 0;                 % Initial temperature in all nodes 
T_right  = -20;        % temperature at x=Lx
T_left   = 20;        % temperature at x=0 

T_top  = 0;             % temperature at z=Lz 
T_bott  = 0;            % temperature at z=0  

T_front=0;            % temperature at y=0 
T_back=0;            % temperature at y=Ly

tol_ss=0.01;            % tolerance for steady state

%% Initial conditions for Steady State
T=zeros(Ny,Nx,Nz);

T(:,1,:)    =  T_left;
T(:,Nx,:)   =  T_right;

T(1,:,:)    =  T_front;
T(Ny,:,:)   =  T_back;

T(:,:,1)    =  T_bott;
T(:,:,Nz)   =  T_top;

%% Heat fluxes at top and bottom (2nd type -- von Neuman)

% enable calculation of temperature
% from heat fluxes at the walls (2nd type B.C. at the walls)
BC2=1; % set to 0 to disable

if BC2==1 % see also lines 130++
    Q_top   =  0; % Q=0 means thermal insulation
    Q_bott  =  0;

    Q_left  =  0;
    Q_right =  0;

    Q_front =  0;
    Q_back  =  0;
end
%% Constants

iter=1;  % counter of iterations
error  = 1; % initial error for iterations

tmax=3600; % max time, sec 
t=0;
k1=alph*dt; 

%% plotting 
T_min=T_right;
T_max=T_left;

numisosurf=25; % number of isosurfaces
isovalues=linspace(T_min,T_max,numisosurf);

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
z=linspace(0,Lz,Nz);
[X Y Z]=meshgrid(x,y,z);

figure('units','pixels','position',[100 100 1280 720])
 
%% Steady-State Solution
tic
i=2:Ny-1; %inner nodes along y
j=2:Nx-1; %inner nodes along x
k=2:Nz-1; %inner nodes along z
T_new=T;

method='steadystate';

while t<=tmax%error >= tol_ss % by tmax or by tolerance (for Staeady State)
  
%% von Neuman B.C. (fluxes at walls)
% grad T as central difference    
if BC2==1
      
% grad T via central difference -- Neuman B.C. at the walls (constant fluxes or grad(T) )
% uncomment the respective walls to enable T calculation at the boundary nodes from Q

        T(1,:,:)  =  T(3,:,:)     +Q_front*2*dy;     %  flux on front
        T(Ny,:,:) =  T(Ny-2,:,:)  +Q_back*2*dy;      %  flux on back
      
%       T(:,1,:)  =  T(:,3,:)    +Q_left*2*dx;       %  flux on left
%       T(:,Nx,:) =  T(:,Nx-2,:) +Q_right*2*dx;      %  flux on right
    
%       T(:,:,1)  =  T(:,:,3)    +Q_bott*2*dz;       %  flux on bottom
       T(:,:,Nz) =  T(:,:,Nz-2) +Q_top*2*dz;        %  flux on top
end    
   

%%     
  switch method
        case 'steadystate'
    
   % Solving Steady-State Laplace's eqn on 3D FD stencil
    T_new(i,j,k)=(dz^2*(T(i,j,k-1)+T(i,j,k+1))+dy^2*(T(i+1,j,k)+T(i-1,j,k))+dx^2*(T(i,j+1,k)+T(i,j-1,k)))/(2*(dx^2+dy^2+dz^2));    
  
        case 'Euler'
    
    % FTCS transient -- forward Euler
   T_new(i,j,k) =T(i,j,k)+k1*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dy^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dx^2)+((T(i,j,k-1)-2*T(i,j,k)+T(i,j,k+1))/dz^2));
 
  end
    iter=iter+1;
    t=iter*dt; % pointless for Steady State
%     error=max(max(max(abs(T_new-T))));    
  T(i,j,k)=T_new(i,j,k);

if mod(iter,2)==0  % e.g. plot every 5 iter
    plot3D(X,Y,Z,Lx,Ly,Lz,dx,dy,dz,T,T_min,T_max,isovalues,iter,t); 
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
























