
function plot3D(X,Y,Z,Lx,Ly,Lz,dx,dy,dz,T,T_min,T_max,isovalues,iter, t);
clf

numisosurf=25; % number of isosurfaces

%% recalculate scales (disable this section to keep the static scale)
%  T_max = max(T,[],'all');
%  T_min = min(T,[],'all');
%  isovalues=linspace(T_min,T_max,numisosurf);

%% colorbar
format bank
cbh = colorbar ; %Create Colorbar
cbh.Ticks = 1:25; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(round(linspace(T_min,T_max,numisosurf),2));    %Replace the labels of these 8 ticks with the numbers 1 to 8
 
%%
num=numel(isovalues);
 for i=1:num     
     p(i)=patch(isosurface(X,Y,Z,T,isovalues(i)));
     isonormals(X,Y,Z,T,p(i));
     set(p(i), 'CData',i); 
hold on
 end
set(p, 'CDataMapping','direct', 'FaceColor','flat', 'EdgeColor','none')
clr = jet(num);
colormap(clr);
caxis([0 num])
%colorbar

box off; grid on; axis tight; daspect([1 1 1])

camproj perspective
camlight; lighting gouraud; alpha(0.75);
axis([0 Lx 0 Ly 0 Lz]); % plot only half along Y to see inside
axis manual
%view(iter,20); % orbiting camera
%view(180,0); %switch on to have static view
view(3) % iso view
title(['Steady State Solution, iter = ' num2str(iter), ' time = ' num2str(t) ,' s']);
drawnow

