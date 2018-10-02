%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                   CREATE CIRCULAR OR ELLIPTICAL BLOBS                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 2nd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% GEO = PP_CreateBlobs(P, r, dr);                                         %
% GEO = PP_CreateBlobs(P, r, dx, dy);                                     %
% GEO = PP_CreateBlobs(P, a, b, dx, dy);                                  %
% GEO = PP_CreateBlobs(P, a, b, da, db, phi);                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function creates blobs of points which can be used later for       %
% particle advection purposes. The function allows the user to create     %
% blobs through several options. In the first option, the user can create %
% circular blobs with equal spacing in the x and y directions. In the     %
% second option, the user can create circular blobs with different grid   %
% spacings in the x and y directions. In the third option, the user can   %
% create ellipses with selected spacing in the x and y directions. In the %
% fourth and final option, the user can create ellipses as in the third   %
% option with various inclinations.                                       %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'a'          - REAL VECTOR                                              %
%              - List of semi-major axis values in the case of elliptical %
%                blobs.                                                   %
% ----------------------------------------------------------------------- %
% 'b'          - REAL VECTOR                                              %
%              - List of semi-minor axis values in the case of elliptical %
%                blobs.                                                   %
% ----------------------------------------------------------------------- %
% 'da'         - REAL VECTOR                                              %
%              - List of spacings in the semi-major axis direction in the %
%                case of a rotated blob.                                  %
% ----------------------------------------------------------------------- %
% 'db'         - REAL VECTOR                                              %
%              - List of spacings in the semi-minor axis direction in the %
%                case of a rotated blob.                                  %
% ----------------------------------------------------------------------- %
% 'dr'         - REAL VECTOR                                              %
%              - List of spacings in the both directions in the case of a %
%                circular blob.                                           %
% ----------------------------------------------------------------------- %
% 'dx'         - REAL VECTOR                                              %
%              - List of spacings in the x direction in the case of a     %
%                simple circular or elliptical blob.                      %
% ----------------------------------------------------------------------- %
% 'dy'         - REAL VECTOR                                              %
%              - List of spacings in the y direction in the case of a     %
%                simple circular or elliptical blob.                      %
% ----------------------------------------------------------------------- %
% 'GEO'        - STRUCT                                                   %
%              - Struct containing fields X and Y holding the positions   %
%                of all the points within all the blobs.                  %
% ----------------------------------------------------------------------- %
% 'P'          - TWO-COLUMN DOUBLE ARRAY                                  %
%              - List of the centers of the blobs with the first column   %
%                representing the x coordinates and the second column the %
%                y coordinates.                                           %
% ----------------------------------------------------------------------- %
% 'phi'        - REAL VECTOR                                              %
%              - List of rotation angles in radians (from 0 to pi) in the %
%                case of rotated blobs.                                   %
% ----------------------------------------------------------------------- %
% 'r'          - REAL VECTOR                                              %
%              - List of radii values in the case of circular blobs.      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Create a set of elliptical blobs with centers (0.4,0.8), (0.4,0.2),     %
% (1.6,0.8) and (1.6, 0.2), all having semi-major axes of 0.1 and semi-   %
% minor axes of 0.05. Use 0.0005 for the grid spacing in both directions. %
% Also, create the ellipses with the following angles (in the order in    %
% which their centers were stated) 0, pi/4, pi/2 and 3*pi/4. Plot the     %
% ellipse positions with a background double gyre flow. Use the domain    %
% (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 and a time     %
% interval of [0,20] with time-step size 0.1. Use A = 0.1, epsilon =      %
% 0.25, and omega = 2*pi/10.                                              %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> quiver(VEC{1}.X(1:4:end,1:4:end), VEC{1}.Y(1:4:end,1:4:end), ...     %
%           VEC{1}.U(1:4:end,1:4:end), VEC{1}.V(1:4:end,1:4:end), 'k');   %
% >> axis([0 2 0 1]);                                                     %
% >> hold on;                                                             %
% >> grd = PP_CreateBlobs([0.4 0.8; 0.4 0.2; 1.6 0.8; 1.6 0.2], 0.1, ...  %
%                         0.05, 0.0005, 0.0005, [0; pi/4; pi/2; 3*pi/4]); %
% >> scatter(grd.X, grd.Y, 'r', 'filled');                                %
% >> hold off;                                                            %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Advect a set of elliptical blobs in a time-dependent double gyre over   %
% the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01     %
% over the time interval [0,20] with time-step size 0.1. Use A = 0.1,     %
% epsilon = 0.25, and omega = 2*pi/10. Use the time-step dt = 0.1 for the %
% advection (i.e. use the 'singlestep' option). Place the initial blobs   %
% centers at (0.4,0.8), (0.4,0.2), (1.6,0.8) and (1.6, 0.2), all having   %
% semi-major axes of 0.1 and semi-minor axes of 0.05. Use 0.001 for the   %
% grid spacing in both directions. Also, create the ellipses with the     %
% following angles (in the order in which their centers were stated) 0,   %
% pi/4, pi/2 and 3*pi/4.                                                  %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> grd = PP_CreateBlobs([0.4 0.8; 0.4 0.2; 1.6 0.8; 1.6 0.2], 0.1, ...  %
%                         0.05, 0.001, 0.001, [0; pi/4; pi/2; 3*pi/4]);   %
% >> A = PP_AdvectGrid(grd, VEC, 0.1, 'singlestep');                      %
% >> for k = 1:length(t)                                                  %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'k');      %
% axis([0 2 0 1]);                                                        %
% hold on;                                                                %
% scatter(mat2vec(A{k}.X), mat2vec(A{k}.Y), 'r', 'filled');               %
% hold off;                                                               %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% mat2vec                                                                 %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_CreateBlobs
function [GEO] = PP_CreateBlobs(P, a, varargin)

if nargin == 3
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, r = a; else r = a(k); end
        if length(varargin{1}) == 1, dr = varargin{1}; else ...
                                     dr = varargin{1}(k); end
        
        x = (P(k,1)-r):dr:(P(k,1)+r);
        y = (P(k,2)+r):-dr:(P(k,2)-r);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        C = ((TMP.X - P(k,1)).^2 + (TMP.Y - P(k,2)).^2 <= r^2);
        C = mat2vec(C.');
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        TMP.X(~C) = [];
        TMP.Y(~C) = [];
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
        
    end
    
elseif nargin == 4
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, r = a; else r = a(k); end
        if length(varargin{1}) == 1, dx = varargin{1}; else ...
                                     dx = varargin{1}(k); end
        if length(varargin{2}) == 1, dy = varargin{2}; else ...
                                     dy = varargin{2}(k); end
        
        x = (P(k,1)-r):dx:(P(k,1)+r);
        y = (P(k,2)+r):-dy:(P(k,2)-r);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        C = ((TMP.X - P(k,1)).^2 + (TMP.Y - P(k,2)).^2 <= r^2);
        C = mat2vec(C.');
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        TMP.X(~C) = [];
        TMP.Y(~C) = [];
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
    
    end
    
elseif nargin == 5
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, r = a; else r = a(k); end
        if length(varargin{1}) == 1, b = varargin{1}; else ...
                                     b = varargin{1}(k); end
        if length(varargin{2}) == 1, da = varargin{2}; else ...
                                     da = varargin{2}(k); end
        if length(varargin{3}) == 1, db = varargin{3}; else ...
                                     db = varargin{3}(k); end
        
        x = (P(k,1)-r):da:(P(k,1)+r);
        y = (P(k,2)+b):-db:(P(k,2)-b);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        C = (((TMP.X - P(k,1))/r).^2 + ((TMP.Y - P(k,2))/b).^2 <= 1);
        C = mat2vec(C.');
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        TMP.X(~C) = [];
        TMP.Y(~C) = [];
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
        
    end
    
elseif nargin == 6
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, r = a; else r = a(k); end
        if length(varargin{1}) == 1, b = varargin{1}; else ...
                                     b = varargin{1}(k); end
        if length(varargin{2}) == 1, da = varargin{2}; else ...
                                     da = varargin{2}(k); end
        if length(varargin{3}) == 1, db = varargin{3}; else ...
                                     db = varargin{3}(k); end
        if length(varargin{4}) == 1, phi = varargin{4}; else ...
                                     phi = varargin{4}(k); end
        
        x = (P(k,1)-r):da:(P(k,1)+r);
        y = (P(k,2)+b):-db:(P(k,2)-b);
        [TMP1.X, TMP1.Y] = meshgrid(x,y);
        TMP.X = (TMP1.X - P(k,1))*cos(phi) ...
              + (TMP1.Y - P(k,2))*sin(phi) ...
              + P(k,1);
        TMP.Y = (TMP1.X - P(k,1))*sin(phi) ...
              - (TMP1.Y - P(k,2))*cos(phi) ...
              + P(k,2);
        C = ((((TMP1.X - P(k,1))*cos(phi) ...
          +    (TMP1.Y - P(k,2))*sin(phi))/r).^2 ...
          +  (((TMP1.X - P(k,1))*sin(phi) ...
          -    (TMP1.Y - P(k,2))*cos(phi))/b).^2 <= 1);
        C = mat2vec(C.');
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        TMP.X(~C) = [];
        TMP.Y(~C) = [];
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
        
    end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A
% Reason(s)
% * N/A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
