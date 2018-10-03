%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                   CREATE SQUARE OR RECTANGULAR PATCHES                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 3rd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2018 Giuseppe Di Labbio                                   %
%                                                                         %
% This program is free software: you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation, either version 3 of the License, or (at your  %
% option) any later version.                                              %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        %
% General Public License for more details.                                %
%                                                                         %
% You should have received a copy of the GNU General Public License along %
% with this program. If not, see <https://www.gnu.org/licenses/>.         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% GEO = PP_CreatePatches(P, s, ds);                                       %
% GEO = PP_CreatePatches(P, s, dx, dy);                                   %
% GEO = PP_CreatePatches(P, a, b, dx, dy);                                %
% GEO = PP_CreatePatches(P, a, b, da, db, phi);                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function creates patches of points which can be used later for     %
% particle advection purposes. The function allows the user to create     %
% patches through several options. In the first option, the user can      %
% create square patches with equal spacing in the x and y directions. In  %
% the second option, the user can create square patches with different    %
% grid spacings in the x and y directions. In the third option, the user  %
% can create rectangles with selected spacing in the x and y directions.  %
% In the fourth and final option, the user can create rectangles as in    %
% the third option with various inclinations.                             %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'a'          - REAL VECTOR                                              %
%              - List of semi-major side length values in the case of     %
%                rectangular patches.                                     %
% ----------------------------------------------------------------------- %
% 'b'          - REAL VECTOR                                              %
%              - List of semi-minor side length values in the case of     %
%                rectangular patches.                                     %
% ----------------------------------------------------------------------- %
% 'da'         - REAL VECTOR                                              %
%              - List of spacings in the semi-major side direction in the %
%                case of a rotated patch.                                 %
% ----------------------------------------------------------------------- %
% 'db'         - REAL VECTOR                                              %
%              - List of spacings in the semi-minor side direction in the %
%                case of a rotated patch.                                 %
% ----------------------------------------------------------------------- %
% 'ds'         - REAL VECTOR                                              %
%              - List of spacings in the both directions in the case of a %
%                square patch.                                            %
% ----------------------------------------------------------------------- %
% 'dx'         - REAL VECTOR                                              %
%              - List of spacings in the x direction in the case of a     %
%                simple square or rectangular patch.                      %
% ----------------------------------------------------------------------- %
% 'dy'         - REAL VECTOR                                              %
%              - List of spacings in the y direction in the case of a     %
%                simple square or rectangular patch.                      %
% ----------------------------------------------------------------------- %
% 'GEO'        - STRUCT                                                   %
%              - Struct containing fields X and Y holding the positions   %
%                of all the points within all the patches.                %
% ----------------------------------------------------------------------- %
% 'P'          - TWO-COLUMN DOUBLE ARRAY                                  %
%              - List of the centers of the patches with the first column %
%                representing the x coordinates and the second column the %
%                y coordinates.                                           %
% ----------------------------------------------------------------------- %
% 'phi'        - REAL VECTOR                                              %
%              - List of rotation angles in radians (from 0 to pi) in the %
%                case of rotated patches.                                 %
% ----------------------------------------------------------------------- %
% 's'          - REAL VECTOR                                              %
%              - List of side length values in the case of square         %
%                patches.                                                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Create a set of rectangular patches with centers (0.4,0.8), (0.4,0.2),  %
% (1.6,0.8) and (1.6, 0.2), all having semi-major side lengths of 0.1 and %
% semi-minor side lengths of 0.05. Use 0.0005 for the grid spacing in     %
% both directions. Also, create the rectangles with the following angles  %
% (in the order in which their centers were stated) 0, pi/4, pi/2 and     %
% 3*pi/4. Plot the rectangle positions with a background double gyre      %
% flow. Use the domain (x,y) = [0,2]x[0,1] with a constant grid spacing   %
% of 0.01 and a time interval of [0,20] with time-step size 0.1. Use A =  %
% 0.1, epsilon = 0.25, and omega = 2*pi/10.                               %
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
% >> grd = PP_CreatePatches([0.4 0.8; 0.4 0.2; 1.6 0.8; 1.6 0.2], ...     %
%                           0.1, 0.05, 0.0005, 0.0005, ...                %
%                           [0; pi/4; pi/2; 3*pi/4]);                     %
% >> scatter(grd.X, grd.Y, 'r', 'filled');                                %
% >> hold off;                                                            %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Advect a set of rectangular patches in a time-dependent double gyre     %
% over the domain (x,y) = [0,2]x[0,1] with a constant grid spacing of     %
% 0.01 over the time interval [0,20] with time-step size 0.1. Use A =     %
% 0.1, epsilon = 0.25, and omega = 2*pi/10. Use the time-step dt = 0.1    %
% for the advection (i.e. use the 'singlestep' option). Place the initial %
% patch  centers at (0.4,0.8), (0.4,0.2), (1.6,0.8) and (1.6, 0.2), all   %
% having semi-major side lengths of 0.1 and semi-minor side lengths of    %
% 0.05. Use 0.001 for the grid spacing in both directions. Also, create   %
% the rectangular pathes with the following angles (in the order in which %
% their centers were stated) 0, pi/4, pi/2 and 3*pi/4.                    %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A    = 0.1;                                                          %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> grd = PP_CreatePatches([0.4 0.8; 0.4 0.2; 1.6 0.8; 1.6 0.2], ...     %
%                           0.1, 0.05, 0.001, 0.001, ...                  %
%                           [0; pi/4; pi/2; 3*pi/4]);                     %
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

%% PP_CreatePatches
function [GEO] = PP_CreatePatches(P, a, varargin)

if nargin == 3
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, s = a; else s = a(k); end
        if length(varargin{1}) == 1, ds = varargin{1}; else ...
                                     ds = varargin{1}(k); end
        
        x = (P(k,1)-s):ds:(P(k,1)+s);
        y = (P(k,2)+s):-ds:(P(k,2)-s);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
        
    end
    
elseif nargin == 4
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, s = a; else s = a(k); end
        if length(varargin{1}) == 1, dx = varargin{1}; else ...
                                     dx = varargin{1}(k); end
        if length(varargin{2}) == 1, dy = varargin{2}; else ...
                                     dy = varargin{2}(k); end
        
        x = (P(k,1)-s):dx:(P(k,1)+s);
        y = (P(k,2)+s):-dy:(P(k,2)-s);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
    
    end
    
elseif nargin == 5
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, s = a; else s = a(k); end
        if length(varargin{1}) == 1, b = varargin{1}; else ...
                                     b = varargin{1}(k); end
        if length(varargin{2}) == 1, da = varargin{2}; else ...
                                     da = varargin{2}(k); end
        if length(varargin{3}) == 1, db = varargin{3}; else ...
                                     db = varargin{3}(k); end
        
        x = (P(k,1)-s):da:(P(k,1)+s);
        y = (P(k,2)+b):-db:(P(k,2)-b);
        [TMP.X, TMP.Y] = meshgrid(x,y);
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        
        GEO.X = [GEO.X; TMP.X];
        GEO.Y = [GEO.Y; TMP.Y];
        
    end
    
elseif nargin == 6
    
    GEO = struct('X', [], 'Y', []);
    for k = 1:size(P,1)
        
        if length(a) == 1, s = a; else s = a(k); end
        if length(varargin{1}) == 1, b = varargin{1}; else ...
                                     b = varargin{1}(k); end
        if length(varargin{2}) == 1, da = varargin{2}; else ...
                                     da = varargin{2}(k); end
        if length(varargin{3}) == 1, db = varargin{3}; else ...
                                     db = varargin{3}(k); end
        if length(varargin{4}) == 1, phi = varargin{4}; else ...
                                     phi = varargin{4}(k); end
        
        x = (P(k,1)-s):da:(P(k,1)+s);
        y = (P(k,2)+b):-db:(P(k,2)-b);
        [TMP1.X, TMP1.Y] = meshgrid(x,y);
        TMP.X = (TMP1.X - P(k,1))*cos(phi) ...
              + (TMP1.Y - P(k,2))*sin(phi) ...
              + P(k,1);
        TMP.Y = (TMP1.X - P(k,1))*sin(phi) ...
              - (TMP1.Y - P(k,2))*cos(phi) ...
              + P(k,2);
        TMP.X = mat2vec(TMP.X.');
        TMP.Y = mat2vec(TMP.Y.');
        
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
