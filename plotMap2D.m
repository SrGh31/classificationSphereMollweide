function [plotPts]= plotMap2D(pseudo_prots,allData,X,Y,Z,projType)    
    opt=projType; % option 3= Mollweide projection
    [ppxn,ppyn,ppzn]=perspective_change(pseudo_prots(:,1),pseudo_prots(:,2),pseudo_prots(:,3),X,Y,Z);
    [azimuth_pp,elevation_pp,rad_pp] = cart2sph(ppxn,ppyn,ppzn);
    [X_pp,Y_pp]=map_projection(opt,azimuth_pp, elevation_pp,rad_pp);
    
    [xn,yn,zn]=perspective_change(allData(:,1),allData(:,2),allData(:,3),X,Y,Z);
    [azimuth_r,elevation_r,rad] = cart2sph(xn,yn,zn);
    [X_r,Y_r]=map_projection(opt,azimuth_r,elevation_r,rad);
    plotPts.real2D=[X_r,Y_r];
    plotPts.prots2D=[X_pp,Y_pp];
    plotPts.real3D=[xn,yn,zn];
    plotPts.prots3D=[ppxn,ppyn,ppzn];
end