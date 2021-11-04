function [X,Y,str] = map_projection(projection, azimuth, elevation,radius, varargin)
% This function asks for the type of projection user wants in argument
% 'projection' and performs that particular projection. The other 2 inputs
% 'azimuth and 'elevation' are the spherical coordinates of the actual 3D
% sphere. projections included in this function for now are: Mercator (1),
% Azimuth:Lambert area-preserving (2), Mollweide (3), and Hammer (4) % projections.

p=inputParser;
p.addRequired('projection',@isnumeric);
p.addRequired('azimuth',@isfloat);
p.addRequired('elevation',@(x)length(x)==length(azimuth));
%p.addOptional('radius',ones(length(elevation),1), @(x) length(x)==length(azimuth));
p.addOptional('azimuth_offset',0, @(x) length(x)==length(azimuth) & @isfloat);
p.addOptional('elevation_offset',0, @(x) length(x)==length(elevation) & @isfloat);
p.parse(projection,azimuth,elevation, varargin{:});


 switch projection
    case 1 %Mercator
        X=radius.*azimuth;
        Y=radius.*(log(tan((pi/4)+0.5*elevation)));
        str='Mercator projection';
    case 2 % Lambert equal area
        X=azimuth+pi/20;
        Y=sin(elevation);
        str='Lambert equal area projection';
    case 3 % Mollweide
        theta=elevation; theta_new=0; count=1; diff=0;
        while (diff>0.00005 || count==1)
            theta_new=theta-bsxfun(@rdivide, (2*theta+sin(2*theta)-pi*sin(elevation)),(2+2*cos(2*theta)));
            diff=mean(abs(theta-theta_new));
            theta=theta_new; 
            count=count+1;
        end
        X=((2^1.5)/pi)*radius.*(azimuth+0.0*ones(length(azimuth),1)).*cos(theta);
        Y=(2^0.5)*radius.*sin(theta);
        str='Mollweide projection';
    case 4 % Hammer
        X=bsxfun(@rdivide, ((2^1.5)*cos(elevation).*sin(azimuth/2)), ((1+(cos(elevation).*cos(azimuth/2)).^0.5)));
        Y=bsxfun(@rdivide, ((2^0.5)*sin(elevation)), ((1+(cos(elevation).*cos(azimuth/2)).^0.5)));
        str='Hammer projection';
 end
end
