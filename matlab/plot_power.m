% serp1_det0;
%sealer_mod4_det0;
%serp1_det0_number2;
smallcrit_det0;

I = DETFApower(:,9:10);

%[X,Y] = meshgrid(DETFApowerCOORD(:,1), DETFApowerCOORD(:,2));

%surf(X, Y, DETFApower(:,11));
 
%plot3(DETFApowerCOORD(:,1), DETFApowerCOORD(:,2),
%DETFApower(:,11), '.');

M = [];
for i = 1:size(DETFApowerCOORD,1)
    if (DETFApower(i,11) > 0.1)
    M = [M DETFApower(i,11)];
    end
end

fprintf('M = %d\n', mean(M));

figure(5); clf
for i = 1:size(DETFApowerCOORD,1)
    %for i = 65:161
  x = DETFApowerCOORD(i,1); 
  y = DETFApowerCOORD(i,2);
  
  %r = 9.155;
  r = 9;
  
  xp = [];
  yp = [];
  
  xp = [xp, x - r];
  xp = [xp, x - r/2];
  xp = [xp, x + r/2];
  xp = [xp, x + r];
  xp = [xp, x + r/2];
  xp = [xp, x - r/2];
  
  yp = [yp, y];
  yp = [yp, y + sqrt(3)/2*r];
  yp = [yp, y + sqrt(3)/2*r];
  yp = [yp, y];
  yp = [yp, y - sqrt(3)/2*r];
  yp = [yp, y - sqrt(3)/2*r];
  
  %xp = x + r*cos(0:5*pi/3);
  %yp = y + r*sin(0:5*pi/3);
  z = DETFApower(i,11)*ones(size(xp))/mean(M);
  
  % if (z > 0.2) | ((x == 0.0) & (y == 0.0))
  %if (z > 0.2) | ((abs(x) < 70.0) & (abs(y) < 70.0))
  xp = xp';
  yp = yp';
  z  = z';
  
  x = xp;
  y = yp;
  
  dt = delaunayTriangulation(x,y) ;
  tri = dt.ConnectivityList ;
  xi = dt.Points(:,1) ; 
  yi = dt.Points(:,2) ; 
  F = scatteredInterpolant(x,y,z);
  zi = F(xi,yi) ;
  trisurf(tri,xi,yi,zi, 'EdgeColor', 'none') 
  hold on
  %end
end
set(gca,'Color','k')
view(2)
shading interp
colormap jet
colorbar
grid off

for i = 1:size(DETFApowerCOORD,1)
%for i = 65:161
  x = DETFApowerCOORD(i,1); 
  y = DETFApowerCOORD(i,2);
  z = DETFApower(i,11)/mean(M);
  
  if (z > 0.2) | ((abs(x) < 70.0) & (abs(y) < 70.0))
  
  text(x,y,2.0, num2str(z,3), 'Color', 'w', 'HorizontalAlignment', ...
       'center', 'FontWeight', 'bold')
  end
end

xlim([-75, 75]);
ylim([-75, 75]);

xlabel('cm', 'interpreter', 'LaTeX')
ylabel('cm', 'interpreter', 'LaTeX')

pbaspect([1 1 1])

figure(1); clf
x = DETFApowerCOORD(:,1); 
y = DETFApowerCOORD(:,2); 
z = DETFApower(:,11)/1.9E6;

dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ; 
yi = dt.Points(:,2) ; 
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi) 
view(2)
shading interp 
colormap jet
colorbar

xlim([-75, 75]);
ylim([-75, 75]);

xlabel('cm', 'interpreter', 'LaTeX')
ylabel('cm', 'interpreter', 'LaTeX')

pbaspect([1 1 1])

%figure(3); clf
%X = reshape(x,[],2);
%Y = reshape(y,[],2);
%Z = reshape(z,[],2);
%surf(X,Y,Z)

figure(2); clf
clear all  
% generate synthetic data
data=[randn(500,2);
    randn(500,1)+3.5, randn(500,1);];
% call the routine, which has been saved in the current directory 
[bandwidth,density,X,Y]=kde2d(data);
% plot the data and the density estimate
contourf(X,Y,density,256, 'LineStyle', 'none'), hold on
%plot(data(:,1),data(:,2),'r.','MarkerSize',5)