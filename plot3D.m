
function plot3D(X,Y,Z,Lx,Ly,Lz,dx,dy,dz,T,isovalues,iter, t);
clf

numisosurf=25; % number of isosurfaces

%% recalculate scales (dsable this section to keep the static scale)
 T_max = max(T,[],'all');
 T_min = min(T,[],'all');
 isovalues=linspace(T_min,T_max,numisosurf);

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
colorbar

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

