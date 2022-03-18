%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            ADVECTION TOOLBOX                            %
%                                                                         %
% advectODEfun                                                            %
% Particle Advection                                                      %
% Advect particles using a MATLAB ODE solver                              %
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
% adv = advectODEfun(coord, vel, tspan, pts);                             %
% adv = advectODEfun(___, Name, Value);                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Advect a number of particles in a flow field using one of MATLAB's      %
% ordinary differential equation solvers.                                 %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2021b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% advectionODE                                                            %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'coord'        STRUCT ARRAY (1 X 1)                                     %
%              ~ Coordinates. The structure array may contain the         %
%                following fields:                                        %
%                1) coord.x representing the space coordinate x;          %
%                2) coord.y representing the space coordinate y;          %
%                3) coord.z representing the space coordinate z;          %
%                4) coord.t representing the time coordinate t.           %
%                Each field is a column vector (one-dimensional array)    %
%                that strictly increases/decreases monotonically from the %
%                first row to the last. Not all fields need be present in %
%                the coordinate structure, however the order must follow  %
%                [x y z t]. For example, the coordinate structure may     %
%                contain fields ordered as [x z t] but not as [t x z].    %
%                Both empty and scalar arrays are permitted.              %
% ----------------------------------------------------------------------- %
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array may contain the      %
%                following fields:                                        %
%                1) vel.u representing the x-component velocity field;    %
%                2) vel.v representing the y-component velocity field;    %
%                3) vel.w representing the z-component velocity field.    %
%                Each field is a two- to four-dimensional array, with the %
%                dimensions corresponding to the order [x y z t]. Not all %
%                fields need be present in the velocity field structure,  %
%                however the order must follow [u v w]. For example, the  %
%                velocity field structure may contain fields ordered as   %
%                [u w] but not as [w u]. Both empty and scalar arrays are %
%                permitted.                                               %
% ----------------------------------------------------------------------- %
% 'tspan'        ARRAY                                                    %
%              ~ Interval of integration. The array can be either a two-  %
%                element vector or a list of desired advection times. In  %
%                the case of a two-element vector [ti tf], the MATLAB ODE %
%                solver returns the solution evaluated at each internal   %
%                integration step within the interval. If more elements   %
%                are specified, the solver returns the particle positions %
%                evaluated only at the given times.                       %
% ----------------------------------------------------------------------- %
% 'pts'          STRUCT (1 X 1)                                           %
%              ~ Initial particle positions. The structure may contain    %
%                the following fields:                                    %
%                1) pts.x representing the x-coordinates of particles;    %
%                2) pts.y representing the y-coordinates of particles;    %
%                3) pts.z representing the z-coordinates of particles.    %
%                Each field is an array of positions. No restriction is   %
%                placed on the construction of the arrays. However, it is %
%                common to either specify a simple list of points or an   %
%                array with same number of spatial dimensions as the      %
%                velocity field.                                          %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'format'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'list'                                          %
%              ~ Specify the output format for the advected particles.    %
%                The option 'list' will organize the output (adv) so that %
%                each column represents a particle pathline. The option   %
%                'grid' will organize the output to match the format of   %
%                the input (pts) with a new dimension representing time.  %
% ----------------------------------------------------------------------- %
% 'interp'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'linear'                                        %
%              ~ Specify the method of interpolation for intermediate     %
%                grid points. See the MATLAB documentation of interpn for %
%                more information.                                        %
% ----------------------------------------------------------------------- %
% 'extrapval'    SCALAR                                                   %
%                Default: NaN                                             %
%              ~ Scalar value assigned to query points falling outside    %
%                the domain. Unlike MATLAB's default for interpn, where   %
%                extrapolation is performed for the 'makima' and 'spline' %
%                methods and NaN is assigned for others, here only scalar %
%                assignment is permitted, regardless of the interpolation %
%                method.                                                  %
% ----------------------------------------------------------------------- %
% 'odeFun'       FUNCTION HANDLE                                          %
%                Default: @ode45                                          %
%              ~ MATLAB ODE solver. The function handle must correspond   %
%                to one of MATLAB's built-in ODE solvers. For example,    %
%                this could be @ode45, @ode23, @ode113, etc.              %
% ----------------------------------------------------------------------- %
% 'odeOpts'      STRUCT (1 X 1)                                           %
%                Default: MATLAB default for specified ODE solver         %
%              ~ Options structure for the selected MATLAB ODE solver.    %
%                See the MATLAB documentation for the specified ODE       %
%                solver for more information.                             %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'adv'          STRUCT (1 X 1)                                           %
%              ~ Advected particle positions. The structure will contain  %
%                the following fields to match the input (pts) as well as %
%                a time field:                                            %
%                1) adv.x representing the x-coordinates of particles;    %
%                2) adv.y representing the y-coordinates of particles;    %
%                3) adv.z representing the z-coordinates of particles;    %
%                4) adv.t representing the t-coordinates of particles.    %
%                Each field is an array of positions. Formatted according %
%                to this 'format' name-value pair. By default, the arrays %
%                are organized such that each column represents a single  %
%                particle's pathline. Alternatively, the arrays can be    %
%                organized to match the format of the input (pts) with an %
%                additional dimension representing time.                  %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Generate an irrotational vortex on the domain (x,y) = ([-1,1],[-1,1])   %
% with a constant grid spacing of 0.01. Use a circulation of Gamma = 5.   %
% Advect a set of 10 particles in time from 0 to 1 time units. Compare    %
% the computed pathlines to the theoretical streamlines using the exact   %
% streamfunction of the irrotational vortex.                              %
%                                                                         %
% >> coord.x       = linspace(-1, 1, 201).';                              %
% >> coord.y       = linspace(-1, 1, 201).';                              %
% >> Gamma         = 5;                                                   %
% >> [vel, ~, cvp] = irrotVortex(coord, Gamma);                           %
% >> coord.t       = linspace(0, 5, 51).';                                %
% >> vel.u         = repmat(vel.u, [1 1 51]);                             %
% >> vel.v         = repmat(vel.v, [1 1 51]);                             %
% >> pts.x         = (0.1:0.1:1).';                                       %
% >> pts.y         = zeros(length(pts.x), 1);                             %
% >> adv           = advectODEfun(coord, vel, [0 1], pts);                %
% >> isovals       = interpn(coord.x, coord.y, imag(cvp), pts.x, pts.y);  %
% >> quiver(coord.x, coord.y, vel.u(:,:,1).', vel.v(:,:,1).', 5, 'k');    %
% >> axis equal tight;                                                    %
% >> hold on;                                                             %
% >> contour(coord.x, coord.y, imag(cvp).', isovals, 'k');                %
% >> plot(adv.x, adv.y, '--r');                                           %
% >> hold off;                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [adv] = advectODEfun(coord, vel, tspan, pts, varargin)


%% PARSE INPUTS

% Input defaults.
default.format    = 'list';
default.interp    = 'linear';
default.extrapval = NaN;
default.odeFun    = @ode45;
default.odeOpts   = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Input checks.
check.tspan     = @(x) validateattributes(x,                            ...
                       {'logical', 'numeric'},                          ...
                       {'finite', 'nonnan', 'real', 'vector'});
check.format    = @(x) any(validatestring(x,                            ...
                       {'list' 'grid'}));
check.interp    = @(x) any(validatestring(x,                            ...
                       {'linear' 'nearest' 'pchip' 'cubic' 'makima'     ...
                        'spline'}));
check.extrapval = @(x) validateattributes(x,                            ...
                       {'logical', 'numeric'},                          ...
                       {'scalar'});
check.odeFun    = @(x) validateattributes(x,                            ...
                       {'function_handle'},                             ...
                       {});
check.odeOpts   = @(x) validateattributes(x,                            ...
                       {'struct'},                                      ...
                       {'size', [1 1]});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord'                                         );
addRequired ( hParser, 'vel'                                           );
addRequired ( hParser, 'tspan'                       , check.tspan     );
addRequired ( hParser, 'pts'                                           );
addParameter( hParser, 'format'   , default.format   , check.format    );
addParameter( hParser, 'interp'   , default.interp   , check.interp    );
addParameter( hParser, 'extrapval', default.extrapval, check.extrapval );
addParameter( hParser, 'odeFun'   , default.odeFun   , check.odeFun    );
addParameter( hParser, 'odeOpts'  , default.odeOpts  , check.odeOpts   );
parse(hParser, coord, vel, tspan, pts, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% ADVECT PARTICLES

% Get field information from the advection initial condition.
fields    = fieldnames(pts);
numFields = length(fields);
szPts     = size(pts.(fields{1}));
Npts      = numel(pts.(fields{1}));

% Reshape the particle position array for advection. Note that to solve for
% multiple particle positions using MATLAB's ODE solvers, the input must be
% formatted as a two-dimensional array where each column represents the set
% of coordinates for one particle.
pts4ode = zeros(numFields, Npts);
for k = 1:numFields
    pts4ode(k,:) = pts.(fields{k})(:).';
end

% Define the advection ordinary differential equation.
ODEdef = @(t,x) advectionODE(t, x, coord, vel, hParser.Results.interp,  ...
                             hParser.Results.extrapval);

% Initialize a structure object to hold the advected particle positions.
for k = 1:numFields, adv.(fields{k}) = 0; end
adv.t = 0;

% Perform the advection.
[adv.t, tmp] = hParser.Results.odeFun(ODEdef, tspan, pts4ode,           ...
                                      hParser.Results.odeOpts);
clear ODEdef pts4ode tspan;

% Reshape the data.
nt = length(adv.t);
switch lower(hParser.Results.format)
    case 'list'
        for k = 1:numFields
            adv.(fields{k}) = tmp(:,k:numFields:end);
        end
    case 'grid'
        for k = 1:numFields
            adv.(fields{k}) = tmp(:,k:numFields:end).';
            adv.(fields{k}) = reshape(adv.(fields{k}), [szPts nt]);
        end
end
clear nt numFields szPts tmp;


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
% 2022/03/18 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
