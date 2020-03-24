% test_my_arrow3d.m

close all;
clear all;

%% Constants

%% Parameters
Xt = [0; 0; 0]; % tail point of vector
Xh = [1; 1; 1]; % head point of vector

arrowR3( Xt, Xh )
hold on;
arrowR3( Xt, [0; 0; 1] )
hold off;

set(gca, 'FontName', 'Times', 'FontSize', 18)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')

