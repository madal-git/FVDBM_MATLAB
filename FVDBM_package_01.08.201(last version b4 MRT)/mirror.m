function mirror_coor = mirror(xy_ori, xy_mirr_1, xy_mirr_2)

% function mirror_coor = mirror(xy_ori, xy_mirr_1, xy_mirr_2) generates the
% [x;y] coordinates of the mirror image of xy_ori based on the mirror
% defined by xy_mirr_1 and xy_mirr_2

if length(xy_ori)~=2
    error('Check the dimension of the coordinates');
end

if length(xy_mirr_1)~=2
    error('Check the dimension of the coordinates');
end

if length(xy_mirr_2)~=2
    error('Check the dimension of the coordinates');
end

xy_center=norm_joint(xy_mirr_2,(xy_mirr_2-xy_mirr_1),xy_ori);

mirror_coor=2*xy_center-xy_ori;