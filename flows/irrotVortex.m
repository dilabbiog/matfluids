%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         PREDEFINED FLOW TOOLBOX                         %
%                                                                         %
% irrotVortex                                                             %
% Potential Flows                                                         %
% Generate the velocity field of an irrotational vortex                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Département de génie mécanique                                          %
% École de technologie supérieure (ÉTS)                                   %
% Montréal, Québec                                                        %
% Canada                                                                  %
%                                                                         %
% Contributors: Giuseppe Di Labbio                                        %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2022 Giuseppe Di Labbio                                   %
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
% vel = irrotVortex(coord, Gamma, a);                                     %
% [vel, cvp] = irrotVortex(coord, Gamma, a);                              %
% [vel, cvp, vgt] = irrotVortex(coord, Gamma, a);                         %
% [___] = irrotVortex(coord, Gamma, a, ctr);                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Generate the velocity field for a classical irrotational (or potential) %
% vortex with strength (or circulation) Gamma. The coordinates of the     %
% vortex centre can also be specified, the default location being at the  %
% origin. The complex velocity potential and velocity gradient tensor may %
% also be generated. The governing equations can be found in nearly any   %
% fluid dynamics textbook such as [1,2].                                  %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% N/A                                                                     %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% [1] Batchelor, G. K. (2007). An Introduction to Fluid Dynamics.         %
%     Cambridge University Press.                                         %
% [2] Munson, B. R., Okiishi, T. H., Huebsch, W. W., & Rothmayer, A. P.   %
%     (2013). Fundamentals of Fluid Mechanics (7th ed.). John Wiley &     %
%     Sons, Inc.                                                          %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'coord'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Coordinates. The structure array contains the following  %
%                fields:                                                  %
%                1) coord.x representing the space coordinate x;          %
%                2) coord.y representing the space coordinate y.          %
%                Each field is a column vector (one-dimensional array)    %
%                that strictly increases/decreases monotonically from the %
%                first row to the last.                                   %
% ----------------------------------------------------------------------- %
% 'Gamma'        REAL FINITE SCALAR                                       %
%              ~ Vortex circulation. Rotational strength (or circulation) %
%                of the irrotational vortex. A positive value corresponds %
%                to counterclockwise rotation.                            %
% ----------------------------------------------------------------------- %
% 'a'            POSITIVE REAL FINITE SCALAR                              %
%              ~ Radius. Radius at which the stream function is set equal %
%                to zero.                                                 %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'ctr'          REAL FINITE TWO-ELEMENT VECTOR                           %
%                Default: [0 0]                                           %
%              ~ Vortex centre. Location of the centre of the vortex      %
%                relative to the origin. The centre location is specified %
%                as [x y].                                                %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array contains the         %
%                following fields:                                        %
%                1) vel.u representing the velocity component along x;    %
%                2) vel.v representing the velocity component along y.    %
%                Each field is a two dimensional array, in ndgrid format  %
%                (i.e., dimensions represent [x y]).                      %
% ----------------------------------------------------------------------- %
% 'cvp'          COMPLEX 2-DIMENSIONAL ARRAY                              %
%              ~ Complex velocity potential. The complex array contains   %
%                the potential as the real part and the streamfunction as %
%                the imaginary part. The output array is in ndgrid format %
%                (i.e., dimensions represent [x y]).                      %
% ----------------------------------------------------------------------- %
% 'vgt'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity gradient tensor. The structure array contains   %
%                the following fields:                                    %
%                1) vgt.ux representing the derivative of u along x;      %
%                2) vgt.uy representing the derivative of u along y;      %
%                3) vgt.vx representing the derivative of v along x;      %
%                4) vgt.vy representing the derivative of v along y.      %
%                Each field is a two dimensional array, in ndgrid format  %
%                (i.e., dimensions represent [x y]).                      %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Generate an irrotational vortex on the domain (x,y) = ([-1,1],[-1,1])   %
% with a constant grid spacing of 0.01. Use a circulation of Gamma = 5    %
% and a zero stream function radius of 1.                                 %
%                                                                         %
% >> coord.x = linspace(-1, 1, 201).';                                    %
% >> coord.y = linspace(-1, 1, 201).';                                    %
% >> Gamma   = 5;                                                         %
% >> a       = 1;                                                         %
% >> vel     = irrotVortex(coord, Gamma, a);                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vel, varargout] = irrotVortex(coord, Gamma, a, varargin)


%% PARSE INPUTS

% Input defaults.
default.ctr = [0 0];

% Input checks.
check.coord = @(x) examineCoord(x, ["x" "y"], 'all');
check.Gamma = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});
check.a     = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'positive', 'scalar'});
check.ctr   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'vector', 'numel', 2});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord' ,               check.coord );
addRequired ( hParser, 'Gamma' ,               check.Gamma );
addRequired ( hParser, 'a'     ,               check.a     );
addOptional ( hParser, 'ctr'   , default.ctr , check.ctr   );
parse(hParser, coord, Gamma, a, varargin{:});
clear check default;

% Additional verifications.
narginchk(3,4);
nargoutchk(0,3);


%% COMPLEX VELOCITY POTENTIAL

% Transform the coordinates.
zc     = hParser.Results.ctr(1) + 1i*hParser.Results.ctr(2);
z0     = coord.x + 1i*coord.y;
z      = z0 - zc;
[X, Y] = ndgrid(real(z), imag(z));
Z      = X + 1i*Y;
dZ     = ones(size(Z));
clear X Y z z0 zc;

% Define the complex velocity potential.
cvp = -1i*(Gamma/2/pi)*log(Z/a);
if nargout > 1, varargout{1} = cvp; end


%% VELOCITY FIELD

% Initialize the velocity field.
vel = struct('u' , 0, 'v' , 0);

% Compute the velocity field.
dcvp  = -1i*(Gamma/2/pi)*dZ./Z;
vel.u =  real(dcvp);
vel.v = -imag(dcvp);
clear cvp dcvp dZ Z;


%% VELOCITY GRADIENT TENSOR

if nargout > 2
    % Initialize the velocity gradient tensor.
    varargout{2} = struct('ux' , 0, 'uy' , 0, 'vx' , 0, 'vy' , 0);

    % Compute the velocity gradient tensor.
    varargout{2}.ux = -(4*pi/Gamma)*vel.u.*vel.v;
    varargout{2}.uy = -(2*pi/Gamma)*vel.v.^2;
    varargout{2}.vx =  (2*pi/Gamma)*vel.u.^2;
    varargout{2}.vy = -varargout{2}.ux;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  NOTES                                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SUPPRESSED MESSAGES                                                     %
%                                                                         %
% Line(s) N/A                                                             %
% Message(s)                                                              %
% * N/A                                                                   %
% Reason(s)                                                               %
% * N/A                                                                   %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% CHANGE LOG                                                              %
%                                                                         %
% 2022/05/19 -- (GDL) Simplified the code, added radius parameter.        %
% 2022/03/18 -- (GDL) Added narginchk.                                    %
% 2022/03/04 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
