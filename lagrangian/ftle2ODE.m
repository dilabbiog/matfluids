%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           LAGRANGIAN TOOLBOX                            %
%                                                                         %
% ftle2ODE                                                                %
% Finite-Time Lyapunov Exponent                                           %
% 2D FTLE field, advection using a MATLAB ODE solver                      %
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
% ftle = ftle2ODE(coord, vel, tset, DT);                                  %
% ftle = ftle2ODE(coord, vel, tset, DT, mask);                            %
% ftle = ftle2ODE(___, Name, Value);                                      %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Compute the finite-time Lyapunov exponent field (cf. [1-6]) for a two-  %
% dimensional velocity field over a specified interval of time and using  %
% a specified integration time. The function permits both forward and     %
% backward integration in time using positive and negative integration    %
% times, respectively.                                                    %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2021b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% advectionODE                                                            %
% advectODEfun                                                            %
% ddxO2uFD                                                                %
% ddxO2uFDg                                                               %
% idNodeidx                                                               %
%                                                                         %
% Acknowledgments:                                                        %
% N/A                                                                     %
%                                                                         %
% References:                                                             %
% [1] Haller, G. (2000). Finding finite-time invariant manifolds in two-  %
%     dimensional velocity fields. Chaos, 10(1), 99-108.                  %
% [2] Haller, G. (2001). Distinguished material surfaces and coherent     %
%     structures in three-dimensional fluid flows. Physica D, 149, 248-   %
%     277.                                                                %
% [3] Shadden, S. C., Lekien, F. & Marsden, J. E. (2005). Definition and  %
%     properties of Lagrangian coherent structures from finite-time       %
%     Lyapunov exponents in two-dimensional aperiodic flows. Physica D,   %
%     212, 271-304.                                                       %
% [4] Green, M. A., Rowley, C. W., & Haller, G. (2007). Detection of      %
%     Lagrangian coherent structures in three-dimensional turbulence.     %
%     Journal of Fluid Mechanics, 572, 111-120.                           %
% [5] Haller, G. (2011). A variational theory of hyperbolic Lagrangian    %
%     coherent structures. Physica D, 240, 574-598.                       %
% [6] Haller, G. (2015). Lagrangian coherent structures. Annual Review of %
%     Fluid Mechanics, 47, 137-162.                                       %
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
% 'tset'         POSITIVE INTEGER TWO-ELEMENT VECTOR                      %
%              ~ Time indices. The range of times over which to compute   %
%                the FTLE field, specified as a two-element vector        %
%                representing the indices in the time coordinate vecotr   %
%                (coord.t). The indices are specified as [t1 t2], where   %
%                t2 >= t1.                                                %
% ----------------------------------------------------------------------- %
% 'DT'           NONZERO REAL SCALAR                                      %
%              ~ Integration time. The time over which to perform the     %
%                advection for each computed FTLE field. Positive values  %
%                will result in forward time integration and negative     %
%                values in backward time integration.                     %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'mask'         LOGICAL ARRAY                                            %
%              ~ Object mask. The array must have the same size as the    %
%                velocity components. Values of zero (false) are used to  %
%                identify nodes outside the flow or region of interest    %
%                while values of one (true) are used to identify nodes    %
%                within the flow.                                         %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'interp'       CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'linear'                                        %
%              ~ Specify the method of interpolation for intermediate     %
%                grid points. See the MATLAB documentation of interpn for %
%                more information.                                        %
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
% ----------------------------------------------------------------------- %
% 'subNaNs'      CASE-INSENSITIVE CHARACTER ARRAY                         %
%                Default: 'off'                                           %
%              ~ Specify whether to substitute NaN values of the FTLE     %
%                field with their last valid value. This option is useful %
%                in cases with fleeing particles.                         %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'ftle'         REAL ARRAY                                               %
%              ~ Finite-time Lyapunov exponent field. The array will have %
%                the same spatial size as the velocity components whereas %
%                its size in time will be given by diff(tset) + 1.        %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Compute the FTLE field of a double gyre at time t = 0 (coord.t(1))      %
% using a forward integration time of 15. Generate the double gyre on the %
% domain (x,y) = ([0,2],[0,1]) using a constant grid spacing of 0.005     %
% over the time interval [0,20] with time step size 0.1. Use A = 0.1,     %
% epsi = 0.25, and omega = 2*pi/10.                                       %
%                                                                         %
% >> coord.x = linspace(0, 2, 401).';                                     %
% >> coord.y = linspace(0, 1, 201).';                                     %
% >> coord.t = linspace(0, 20, 201).';                                    %
% >> A       = 0.1;                                                       %
% >> epsi    = 0.25;                                                      %
% >> omega   = 2*pi/10;                                                   %
% >> vel     = doubleGyre(coord, A, epsi, omega);                         %
% >> ftle    = ftle2ODE(coord, vel, [1 1], 15);                           %
% >> pcolor(coord.x, coord.y, ftle.');                                    %
% >> shading interp;                                                      %
% >> axis equal tight;                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Compute the FTLE field of a double gyre at time t = 20 (coord.t(end))   %
% using a backward integration time of 15. Generate the double gyre on    %
% the domain (x,y) = ([0,2],[0,1]) using a constant grid spacing of 0.005 %
% over the time interval [0,20] with time step size 0.1. Use A = 0.1,     %
% epsi = 0.25, and omega = 2*pi/10.                                       %
%                                                                         %
% >> coord.x = linspace(0, 2, 401).';                                     %
% >> coord.y = linspace(0, 1, 201).';                                     %
% >> coord.t = linspace(0, 20, 201).';                                    %
% >> A       = 0.1;                                                       %
% >> epsi    = 0.25;                                                      %
% >> omega   = 2*pi/10;                                                   %
% >> vel     = doubleGyre(coord, A, epsi, omega);                         %
% >> ftle    = ftle2ODE(coord, vel, [201 201], -15);                      %
% >> pcolor(coord.x, coord.y, ftle.');                                    %
% >> shading interp;                                                      %
% >> axis equal tight;                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ftle] = ftle2ODE(coord, vel, tset, DT, varargin)


%% PARSE INPUTS

% Input defaults.
default.mask    = [];
default.subNaNs = 'off';
default.interp  = 'linear';
default.odeFun  = @ode45;
default.odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Input checks.
check.coord   = @(x) isfield(x, 't');
check.tset    = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'positive', 'real', 'nondecreasing',    ...
                      'vector', 'numel', 2});
check.DT      = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'nonzero', 'real'});
check.mask    = @(x) validateattributes(x,                              ...
                     {'logical'});
check.subNaNs = @(x) any(validatestring(x,                              ...
                     {'on' 'off'}));
check.interp  = @(x) any(validatestring(x,                              ...
                     {'linear' 'nearest' 'pchip' 'cubic' 'makima'       ...
                      'spline'}));
check.odeFun  = @(x) validateattributes(x,                              ...
                     {'function_handle'},                               ...
                     {});
check.odeOpts = @(x) validateattributes(x,                              ...
                     {'struct'},                                        ...
                     {'size', [1 1]});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'coord'                                     );
addRequired ( hParser, 'vel'                                       );
addRequired ( hParser, 'tset'                      , check.tset    );
addRequired ( hParser, 'DT'                        , check.DT      );
addOptional ( hParser, 'mask'    , default.mask    , check.mask    );
addParameter( hParser, 'subNaNs' , default.subNaNs , check.subNaNs );
addParameter( hParser, 'interp'  , default.interp  , check.interp  );
addParameter( hParser, 'odeFun'  , default.odeFun  , check.odeFun  );
addParameter( hParser, 'odeOpts' , default.odeOpts , check.odeOpts );
parse(hParser, coord, vel, tset, DT, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% INITIALIZATIONS

% Determine the size of the grid in space and time.
N.x = length(coord.x);
N.y = length(coord.y);
N.t = length(coord.t);

% Determine the discretization size in space and time.
d.x = diff(coord.x(1:2));
d.y = diff(coord.y(1:2));
d.t = diff(coord.t(1:2));

% Initialize the FTLE array.
ftle = zeros(N.x, N.y, diff(tset)+1);

% Generate a grid of particles.
[pts.x, pts.y] = ndgrid(coord.x, coord.y);

% Adjust variables for backward time integration.
tvec = (tset(1):tset(2)).';
if (DT < 0)
    coord.t = flipud(coord.t);
    vel.u   = flip(vel.u, 3);
    vel.v   = flip(vel.v, 3);
    tvec    = (N.t + 1) - tvec;
end

if strcmpi(hParser.Results.subNaNs, 'on'), tspan = DT*[0 1];
else,                                      tspan = DT*[0 0.5 1];
end

idx = 0;
for kt = tvec
    idx = idx + 1;
    
    % Perform the advection.
    tadv = coord.t(kt) + tspan;
    if isempty(hParser.Results.mask)
        adv = advectODEfun(coord, vel, tadv, pts,                       ...
                          'format' , 'grid',                            ...
                          'interp' , hParser.Results.interp,            ...
                          'odeFun' , hParser.Results.odeFun,            ...
                          'odeOpts', hParser.Results.odeOpts);
    else
        adv = advectODEfun(coord, vel, tadv, pts,                       ...
                           hParser.Results.mask,                        ...
                          'format' , 'grid',                            ...
                          'interp' , hParser.Results.interp,            ...
                          'odeFun' , hParser.Results.odeFun,            ...
                          'odeOpts', hParser.Results.odeOpts);
    end
    
    if strcmpi(hParser.Results.subNaNs, 'on')
        % Permute the result of the advection so that the first dimension
        % represents time.
        adv.x = permute(adv.x, [3 1 2]);
        adv.y = permute(adv.y, [3 1 2]);
        
        % Determine the indices at which NaN values appear.
        [idxl, ~, idxt, idxs] = idNodeidx(isnan(adv.x));

        % Set all NaN values at the final time to their last valid values.
        if ~isempty(idxl)
            idxl(idxl == 1) = 2;
            adv.x(idxt) = adv.x(idxl-1);
            adv.y(idxt) = adv.y(idxl-1);
        end
        if ~isempty(idxs)
            adv.x(idxs) = adv.x(idxs-1);
            adv.y(idxs) = adv.y(idxs-1);
        end
        
        % Inverse permute the advection result.
        adv.x = ipermute(adv.x, [3 1 2]);
        adv.y = ipermute(adv.y, [3 1 2]);
    end

    % Keep only the final integration time.
    adv.x = adv.x(:,:,end);
    adv.y = adv.y(:,:,end);
    
    if isempty(hParser.Results.mask)
        dgt11 = ddxO2uFD(adv.x(:,:,end), d.x, 1);
        dgt12 = ddxO2uFD(adv.x(:,:,end), d.y, 2);
        dgt21 = ddxO2uFD(adv.y(:,:,end), d.x, 1);
        dgt22 = ddxO2uFD(adv.y(:,:,end), d.y, 2);
    else
        mask  = hParser.Results.mask(:,:,end);
        dgt11 = ddxO2uFDg(adv.x(:,:,end), d.x, mask, 1);
        dgt12 = ddxO2uFDg(adv.x(:,:,end), d.y, mask, 2);
        dgt21 = ddxO2uFDg(adv.y(:,:,end), d.x, mask, 1);
        dgt22 = ddxO2uFDg(adv.y(:,:,end), d.y, mask, 2);
    end
    
    % Compute the elements of the right Cauchy-Green strain tensor.
    cg11 = dgt11.^2 + dgt21.^2;
    cg12 = dgt11.*dgt12 + dgt21.*dgt22;
    cg22 = dgt12.^2 + dgt22.^2;
    
    % Compute elements of the quadratic formula.
    b = cg11 + cg22;
    c = cg11.*cg22 - cg12.^2;
    d = sqrt(b.^2 - 4*c);

    % Compute the FTLE field.
    ftle(:,:,idx) = log(0.5*max(b + d, b - d))/(2*abs(DT));
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
% 2022/06/17 -- (GDL) Simplified the code to use (x,y,t) and (u,v).       %
% 2022/06/17 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% [1] Complete the mask support functionality.                            %
% [2] Find a more clever way to handle different coordinate combinations. %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
