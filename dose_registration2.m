function [cube_d_new,dose_info] = dose_registration2(plan_path,Dose_str)
[~, xV, yV, zV] = load_ct_cube(plan_path);
[dA, xVD, yVD, zVD, dose_info] = load_dose_cube(Dose_str);
[x, y, z] = meshgrid(xV,yV,zV);
delta = 1e-8;
zVD(1) = zVD(1)-1e-3;
zVD(end) = zVD(end)+1e-3;
cube_d_new = finterp3_m(x, y, z, dA, [xVD(1)-delta xVD(2)-xVD(1) xVD(end)+delta], [yVD(1)+delta yVD(2)-yVD(1) yVD(end)-delta], zVD, 0);
cube_d_new(isnan(cube_d_new)) = 0;
cube_d_new(cube_d_new<0) = 0;

% try below if above not working
% [~, xVec, yVec, zVec] = load_ct_cube(plan_path);
% [cube_d, xVec_d, yVec_d, zVec_d, dose_info] = load_dose_cube(plan_path,Dose_str);
% [x, y, z] = meshgrid(yVec_d,xVec_d,zVec_d);
% 
% [xi, yi, zi] = meshgrid(yVec,xVec,zVec);
% cube_d_new = interp3(x,y,z,cube_d,xi,yi,zi,'linear',0);
% 
% cube_d_new(isnan(cube_d_new)) = 0;
% cube_d_new(cube_d_new<0) = 0;
% 
% cube_d_new = rot90(cube_d_new);
% cube_d_new = rot90(cube_d_new);
% cube_d_new = flip(cube_d_new);