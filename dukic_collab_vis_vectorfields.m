%% Xeq view
x = linspace(0,0.01,25);
y = x';
[xx,yy] = ndgrid(x,y);

Y = [0.01 0.01];
Z = arrayfun(@(x,y)rhs_X([x y],Y),xx,yy,'uni',0);
Z_true = arrayfun(@(x,y)rhs_X_true([x y],Y),xx,yy,'uni',0);

subplot(2,1,1)
surf(xx,yy,cellfun(@(z)z(2),Z'),'edgecolor','none')
view([0 90])
colorbar
subplot(2,1,2)
surf(xx,yy,cellfun(@(z)z(2),Z_true'),'edgecolor','none')
view([0 90])
colorbar
%% Yeq view
x = linspace(0.001,0.1,50);
y = x';
[xx,yy] = ndgrid(x,y);

X = [0.1 0.1];
Z = arrayfun(@(x,y)rhs_Y(X,[x y]),xx,yy,'uni',0);
Z_true = arrayfun(@(x,y)rhs_Y_true(X,[x y]),xx,yy,'uni',0);

subplot(2,1,1)
surf(xx,yy,cellfun(@(z)z(1),Z'),'edgecolor','none')
view([0 90])
colorbar
subplot(2,1,2)
surf(xx,yy,cellfun(@(z)z(1),Z_true'),'edgecolor','none')
view([0 90])
colorbar
