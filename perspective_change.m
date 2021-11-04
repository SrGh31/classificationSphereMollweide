function [X_r,Y_r,Z_r] = perspective_change(X,Y,Z,A,B,C)
% perspective changes the viewpoint of the 3D matrix to one specified by
% angles A, B, and C (A,B,C are in degrees)
%Conversion of A,B,C to radians
A=A/180*pi; B=B/180*pi; C=C/180*pi;
dat=[X Y Z];
%Defining the rotaion matrices
R_x=[1         0         0; 0      cos(A) sin(A); 0      -sin(A) cos(A)];
R_y=[cos(B)    0   -sin(B); 0        1         0; sin(B)  0      cos(B)];
R_z=[cos(C) sin(C)       0;-sin(C) cos(C)      0; 0       0           1];
%Performing rotation
for idx=1:length(X)
    dat_x(idx,:)=R_x*dat(idx,:)';
    dat_y(idx,:)=R_y*dat_x(idx,:)';
    dat_z(idx,:)=R_z*dat_y(idx,:)';
end
X_r=dat_z(:,1); Y_r=dat_z(:,2); Z_r=dat_z(:,3);
end

