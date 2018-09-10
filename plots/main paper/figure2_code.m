clc
clear all
x = 0:0.001:4;
y = 0:0.001:4;
[X,Y] = meshgrid(x,y);
Z = sqrt(X.*Y-((X+Y-1)./(2*sqrt(2))).^4);
Z(imag(Z)~=0)=nan;
figure
surfl([X,X],[Y,Y],[Z,-Z])
set(gca, 'FontSize', 12);
xlabel('$S_{xx}$','FontSize', 18, 'Interpreter', 'latex');
ylabel('$S_{yy}$','FontSize', 18, 'Interpreter', 'latex');
zlabel('$S_{xy}$','FontSize', 18, 'Interpreter', 'latex');
colormap(bone)
shading interp
remove_border()
