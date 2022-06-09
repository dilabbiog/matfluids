%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         PREDEFINED FLOW TOOLBOX                         %
%                                                                         %
% uniformFlow                                                             %
% Potential Flows                                                         %
% Generate the velocity field of a uniform flow                           %
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
% vel = uniformFlow(coord, Uinf);                                         %
% [vel, cvp] = uniformFlow(coord, Uinf);                                  %
% [vel, cvp, vgt] = uniformFlow(coord, Uinf);                             %
% [___] = uniformFlow(coord, Uinf, ctr);                                  %
% [___] = uniformFlow(coord, Uinf, ctr, alpha);                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Generate the velocity field for a uniform flow with speed Uinf. A       %
% translation and rotation of the potential flow can also be specified,   %
% the angle being measured counterclockwise from the positive x axis. The %
% default location is at the origin and the default angle is zero. The    %
% complex velocity potential and the velocity gradient tensor may also be %
% generated. The governing equations can be found in nearly any fluid     %
% dynamics textbook such as [1,2].                                        %
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
% 'Uinf'         REAL FINITE SCALAR                                       %
%              ~ Flow speed. Speed of the uniform flow. A positive value  %
%                corresponds to flow in the positive x direction.         %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'ctr'          REAL FINITE TWO-ELEMENT VECTOR                           %
%                Default: [0 0]                                           %
%              ~ Potential flow centre. Location, relative to the origin, %
%                at which the complex velocity potential is zero. The     %
%                centre location is specified as [x y].                   %
% ----------------------------------------------------------------------- %
% 'alpha'        REAL FINITE SCALAR                                       %
%                Default: 0                                               %
%              ~ Potential flow angle. Rotation of the potential flow by  %
%                an angle measured counterclockwise from the positive x   %
%                axis.                                                    %
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
% Generate a uniform flw on the domain (x,y) = ([-1,1],[-1,1]) with a     %
% constant grid spacing of 0.01. Use a flow speed of 1.                   %
%                                                                         %
% >> coord.x = linspace(-1, 1, 201).';                                    %
% >> coord.y = linspace(-1, 1, 201).';                                    %
% >> Uinf    = 1;                                                         %
% >> vel     = uniformFlow(coord, Uinf);                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vel, varargout] = uniformFlow(coord, Uinf, varargin)


%% PARSE INPUTS

% Input defaults.
default.ctr   = [0 0];
default.alpha = 0;

% Input checks.
check.coord = @(x) examineCoord(x, ["x" "y"], 'all');
check.Uinf  = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});
check.ctr   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'vector', 'numel', 2});
check.alpha = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'real', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord' ,                 check.coord );
addRequired ( hParser, 'Uinf'  ,                 check.Uinf  );
addOptional ( hParser, 'ctr'   , default.ctr   , check.ctr   );
addOptional ( hParser, 'alpha' , default.alpha , check.alpha );
parse(hParser, coord, Uinf, varargin{:});
clear check default;

% Additional verifications.
narginchk(2,4);
nargoutchk(0,3);


%% COMPLEX VELOCITY POTENTIAL

% Transform the coordinates.
zc      = hParser.Results.ctr(1) + 1i*hParser.Results.ctr(2);
alpha   = hParser.Results.alpha;
[X0,Y0] = ndgrid(coord.x, coord.y);
Z0      = X0 + 1i*Y0;
Zp      = Z0 - zc;
Z       = Zp*exp(-1i*alpha);
dZ      = exp(-1i*alpha)*ones(size(Z));
clear alpha X0 Y0 Z0 zc;

% Define the complex velocity potential.
cvp = Uinf*Z;
if nargout > 1, varargout{1} = cvp; end

%% VELOCITY FIELD

% Initialize the velocity field.
vel = struct('u' , 0, 'v' , 0);

% Compute the velocity field.
dcvp  = Uinf*dZ;
vel.u =  real(dcvp);
vel.v = -imag(dcvp);
clear cvp dcvp dZ Z;


%% VELOCITY GRADIENT TENSOR

X = real(Zp);
Y = imag(Zp);
if nargout > 2
    % Initialize the velocity gradient tensor.
    varargout{2} = struct('ux' , 0, 'uy' , 0, 'vx' , 0, 'vy' , 0);
end
clear X Y Zp;


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
% 2022/06/09 -- (GDL) Fixed grid error.                                   %
% 2022/06/09 -- (GDL) Added centre variable to set zero potential point.  %
% 2022/05/20 -- (GDL) Style modifications for consistency purposes.       %
% 2022/05/19 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
