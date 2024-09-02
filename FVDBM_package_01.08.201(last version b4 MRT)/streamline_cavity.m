disp('Before running this code, please make sure you have done the following steps:')
disp('1. Run function paraview_stream().Please to refer to the function file to see exact attributes and arguments;')
disp('2. Read Cavity_streamline.mat into the current workspace after clearing the workspace if necessary.')

startx_domain = [0 0.025 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.925 0.95 0.975 1.0];
starty_domain = [0 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
[startx, starty] = meshgrid(startx_domain, starty_domain);
figure;
% streamline(X,Y,u,v,startx,starty,0.01)
streamline(X,Y,U(1,:)',U(2,:)',startx,starty,0.01)
axis equal tight