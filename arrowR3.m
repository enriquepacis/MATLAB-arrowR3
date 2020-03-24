function arrowR3( X_tail, X_head, varargin )
% arrowR3() draws a 3D arrow in a MATLAB plot.
% 
% SYNTAX
%
%  arrowR3( X_tail, X_head ) draws a 3D arrow from point Xi to Xf
%
%  arrowR3( X_tail, X_head, Prop1, Val1, Prop2, Val2, ...) sets additional
%           properties of the 3D arrow through property-value pairs
%
%
% INPUTS
%  Xi, Xf   Xi is the "tail" of the arrow, and Xf is the "head".
%
%
% PROPERTIES
%  The following properties may be modified using property-value pairs.
% 
%  'ArrowHeadRadius'    Radius of the arrowhead at its base. The default
%                       value is 0.05, or 5% of the arrow length.
%
%
%  'ArrowHeadLength'    Length of the arrowhead. The default value is 
%                       0.075, or 7.5% of the arrow length.
%
%  
%  'ShaftRadius'        Radius of the arro shaft. The default value is 
%                       0.0125, or 1.25% of the arrow length.
%
%
%  'Vertices'           Vertices for cylinders used to visualize the arrow.
%                       The default value is 12.
%
%
%  'ArrowHeadColor'     This is an RGB triple of the form [R G B] that
%                       specifies arrowhead color.
%
%                       Default value: [11 131 222]/255
%
%
%  'ArrowHeadBackColor' This is an RGB triple that specifies the color
%                       of the back of the arrowhead. I recommend making
%                       it slightly darker than 'ArrowHeadColor'. This is a
%                       subtle (but very helpful) feature, since it helps
%                       the viewer more clearly see/understand the
%                       direction of the arrow. As a suggestion, one might
%                       set 'ArrowHeadColor' to some [R G B], and then set
%                       'ArrowHeadBackColor' to 0.6*[R G B].
%
%                       Default value: [7 86 145]/255;
%
%
%  'ShaftColor'         This is an RGB triple of the form [R G B] that
%                       specifies arrowhead color.
%
%                       Default value: [145,90,7]/255
%
%
% Version 0.0
% By E.P. Blair
% Baylor University
% Created 2020.03.24 during the coronavirus pandemic
%
% This was inspired by arrow3d() [author: Changshun Deng; dated: August
% 30, 2005]. In MATLAB 2020a, features of Deng's code no longer seem to be
% supported, but this seemed to me to produce the best 3D arrows I've ever
% used.
%
% At this time, arrowR3 supports only one arrow. Todo:
%  - Enable multi-arrow support
%  - Enable support for fixed lengths for arrowhead radius, length, and
%    shaft radius
%

%% Default parameters
relative_lengths = 1; % [bool] toggles whether shaft radius, arrowhead
                      %        radius, and arrowhead length are relative or
                      %        fixed. Presently, fixed lengths are not
                      %        supported.

                      
R_shaft = 0.0125; % cylinder radius relative to overall arrow length
L_ah = 0.075; % arrowhead length relative to overall arrow length
W_ah = 0.05 % width of arrowhead at base, relative to overall arrow length
N = 12; % number of points for surfaces (shaft, arrowhead)
AHeadColor = [11 131 222]/255;
AHeadBackColor = [7 86 145]/255;
ShaftColor = [145,90,7]/255;


%% Default parameter override
args = varargin;
while length(args) >= 2
    prop = args{1};
    val = args{2};
    args = args(3:end);
    
    switch prop
        case 'ShaftRadius'
            R_shaft = val;
        case 'ArrowHeadRadius'
            W_ah = val;
        case 'ArrowHeadLength'
            L_ah = val;
        case 'ArrowHeadColor'
            AHeadColor = val;
        case 'ArrowHeadBackColor'
            AHeadBackColor = val;
        case 'ShaftColor'
            ShaftColor = val;
        otherwise
            error(['ARROWR3: ''', prop, ''' is an invalid property name.'])
    end
    
end

r = X_head - X_tail; % displacement vector from tail to head

[az, el, L] = cart2sph(r(1), r(2), r(3));

% cart2sph returns an elevation referenced from the xy plane, but I want to
% rotate the cylinder objects away from the z axis. The rotation angle,
% therefore, is pi/2 - el, so that if the elevation is 0, then we need a
% pi/2 rotation about the y axis. The azimuthal rotation needs no
% adjustment.
ty = pi/2 - el;
tz = az;
Ry = rotationMatrix('y', ty);
Rz = rotationMatrix('z', tz);
% cylinder returns 2xN matrices of x, y, and z points
%   One row is the bottom face of the cylinder, and the other row is for
%   the top face
%   The bottom face is centered at the origin and lies within the x-y plane
%   The top face is centered at the origin and lies in the plane z = 1.
%   This standard cylinder will be stretched, rotated, and displaced so
%   that the new "bottom" face is at the "tail" point.
[x, y, z] = cylinder(L*R_shaft, N);


% storage matrices for rotated shaft coordinates
xr = zeros(size(x));
yr = zeros(size(y));
zr = zeros(size(z));
for row = 1:2
    Xr = Rz*Ry*[x(row,:); ...
        y(row,:); ...
        L*(1-L_ah)*z(row,:)]; % stretch and rotate
    % the stretching is only to the length L*(1-La_h) to leave room for the
    % arrow head.
    xr(row,:) = Xr(1,:) + X_tail(1); % displace to tail point, Xi
    yr(row,:) = Xr(2,:) + X_tail(2); % displace to tail point, Xi
    zr(row,:) = Xr(3,:) + X_tail(3); % displace to tail point, Xi
end

% Initial standard cylinder for arrowhead, with taper from R = W_ah*L to
% R = 0

[xah, yah, zah] = cylinder(W_ah*L*[1, 0], N);
% storage matrices for arrowhead rotated coordinates
xahr = zeros(size(xah));
yahr = zeros(size(yah));
zahr = zeros(size(zah));

[dx_se, dy_se, dz_se] = sph2cart(az, el, L*(1-L_ah));
shaft_endpoint = X_tail + [dx_se; dy_se; dz_se];
for row = 1:2
    Xahr = Rz*Ry*[xah(row,:); ...
        yah(row,:); ...
        L*L_ah*zah(row,:)]; % stretch and rotate
    
    xahr(row,:) = Xahr(1,:) + shaft_endpoint(1); % displace to head point, Xi
    yahr(row,:) = Xahr(2,:) + shaft_endpoint(2); % displace to head point, Xi
    zahr(row,:) = Xahr(3,:) + shaft_endpoint(3); % displace to head point, Xi
end

hold on;
% surface for the shaft
surf(xr, yr, zr, 'EdgeColor','None','FaceColor', ShaftColor);
% surface for the arrow head
surf(xahr, yahr, zahr, 'EdgeColor','None','FaceColor', AHeadColor)
% patch for the arrow head - back side
patch(xahr(1,:),yahr(1,:),zahr(1,:), AHeadBackColor);
hold off;

axis equal
% lighting gouraud

end

function R = rotationMatrix(varargin)
% rotationMatrix constructs a rotation matrix in dim
% dimensions about an axis
%
% SYNTAX
% rotationMatrix()
%
% rotationMatrix(ax, theta)
%
% rotationMatrix(ax, theta, dim)
% 
% INPUTS
% 
%   dim    [int]     Number of dimensions. {2, [3]}
%
% 'x', 'y', or 'z' of angle theta
%
% By E.P. Blair
% Baylor University
% Created 2020.03.24 during the coronavirus pandemic
%

switch nargin
    case 2
        ax = varargin{1};
        theta = varargin{2};
        dim = 3;
        
    case 3
        ax = varargin{1};
        theta = varargin{2};
        dim = varargin{3};
        
    otherwise
        error('rotationMatrix: invalid number of arguments.')
end
        
        

switch ax
    case 'x'
        R = [1, 0, 0;
            0, cos(theta), -sin(theta);
            0, sin(theta), cos(theta)];
    case 'y'
        R = [cos(theta), 0, sin(theta);
            0, 1, 0;
            -sin(theta), 0, cos(theta)];
    case 'z'
        R = [cos(theta), -sin(theta), 0;
            sin(theta), cos(theta), 0;
            0, 0, 1];
    otherwise
        error(['Invalid axis input: ', ax])
end


end

