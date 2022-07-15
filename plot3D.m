
function plot3D(X,Y,Z,Lx,Ly,Lz,T,isovalues,iter, t);
clf

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
axis([0 Lx Ly/2 Ly 0 Lz]); % plot only half along Y to see inside
axis manual
view(iter*2,20);
%view(3); switch on to have static view
title(['Steady State Solution, iter = ' num2str(iter), ' time = ' num2str(t) ,' s']);
drawnow

