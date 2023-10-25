function Rot_vec_josana_v1(x,y,vx,vy,degree,mask)

% function to rotate the vectors and assign the vectors to the new rotated
% pixel positions

Folder_path = cd;

Rot.degree = degree;
Rot.mfile = 'Rot_vec_josana_v1';

% Change sign of vy so as to bring vx and vy in the right handed coodinate
% system
vy = - (vy);

% get the size of each dimension of velocity matrices
s3 = size(vx,3);
s2 = size(vx,2);
s1 = size(vx,1);

% convert degree to radians
Rot.degree_rad = degree*pi/180;

% rotation matrix to rotate vectors anticlockwise
R = [cos(Rot.degree_rad) -sin(Rot.degree_rad); sin(Rot.degree_rad) cos(Rot.degree_rad)];

% rotation matrix for using in the maketform function
Rot.RM = [cos(Rot.degree_rad) -sin(Rot.degree_rad) 0; sin(Rot.degree_rad) cos(Rot.degree_rad) 0; 0 0 1];
rot = maketform('affine',Rot.RM);

Rot.Rotation_transformation = rot;

[Rot.Irot_thresh,xdata,ydata] = imtransform(mask,rot);
[m n] = find(Rot.Irot_thresh == 1);
Rot.Ppole = max(max(n));
Rot.Apole = min(min(n));
Rot.ELength = Rot.Ppole - Rot.Apole;
Rot.h_top = min(min(m));
Rot.h_bottom = max(max(m));
Rot.ecenter_y = (Rot.h_bottom - Rot.h_top)/2;

% preallocate the final output matrices
Rot.vxrot = nan(s1,s2,s3);
Rot.vyrot = nan(s1,s2,s3);

% Rotate the initial xy coordinates of the velocity vectors
[Rotx,Roty] = tformfwd(rot,x(:,:,1),y(:,:,1));

% Add the new positions to the starting row and column indices to get the
% final true values of the rotated xy coordinates
x_origin = xdata(1) - 0.5;
y_origin = ydata(1) - 0.5;
Rot.xrot = round(Rotx - x_origin);
Rot.yrot = round(Roty - y_origin);

% Rotate the vectors and assign the vectors to the new xy coordinates
for i = 1:s3
    
    for j = 1:s1
        
        for k = 1:s2
            
            Rot_vx_vy = R * [vx(j,k,i); vy(j,k,i)];
            Rot.vxrot(j,k,i) = Rot_vx_vy(1);
            Rot.vyrot(j,k,i) = Rot_vx_vy(2);
        end
    end
end

% clf;subplot(1,2,1),imshow(Rot.Irot_thresh)
% hold on
% quiver(Rot.xrot,Rot.yrot,vx(:,:,10),vy(:,:,10),0,'g')
% subplot(1,2,2),imshow(Rot.Irot_thresh)
% hold on
% quiver(Rot.xrot,Rot.yrot,Rot.vxrot(:,:,10),Rot.vyrot(:,:,10),0,'g')
% saveas(gcf,strcat(Folder_path,'/Analysis_v1/Rotated_quiver'),'tif')

save(strcat(Folder_path,'/Analysis_v1/Rot_pivlab.mat'),'Rot')
% disp('--------Rotation of vectors over---------')


