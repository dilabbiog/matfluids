%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% geospace                                                                %
% Array Generation                                                        %
% Generate an array with geometric series spacing                         %
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
% y = geospace(x1, x2, n);                                                %
% y = geospace(x1, x2, n, 'reverse');                                     %
% y = geospace(___, Name, Value);                                         %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Generate an array (column vector) with a spacing defined by a geometric %
% series between two points. This function complements the linspace and   %
% logspace functions of MATLAB. By default, the spacing decreases from    %
% the first to the last point. By specifying the 'reverse' option, the    %
% spacing will instead increase from the first to the last point. This    %
% function supports variable precision arithmetic.                        %
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
% N/A                                                                     %
%                                                                         %
% ======================================================================= %
% Input Arguments (Required):                                             %
% ----------------------------------------------------------------------- %
% 'x1'           LOGICAL/NUMERIC/VPA REAL SCALAR                          %
%              ~ Array limit 1. The first point in the array that will be %
%                generated.                                               %
% ----------------------------------------------------------------------- %
% 'x2'           LOGICAL/NUMERIC/VPA REAL SCALAR                          %
%              ~ Array limit 2. The last point in the array that will be  %
%                generated.                                               %
% ----------------------------------------------------------------------- %
% 'n'            POSITIVE INTEGER SCALAR                                  %
%              ~ Array size. The number of points to generate the array   %
%                with.                                                    %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'reverse'      CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'reverse' to have the spacing increase from x1   %
%                to x2 rather than decrease (default).                    %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'tol'          NONNEGATIVE REAL SCALAR                                  %
%                Default: 10^-6                                           %
%              ~ Tolerance on the common ratio of the geometric series.   %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'y'            NUMERIC REAL ARRAY                                       %
%              ~ Output array. The generated column vector monotonically  %
%                increasing/decreasing from x1 to x2.                     %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Create an array with spacing decreasing as a geometric series from 2 to %
% 7. Use 51 points.                                                       %
%                                                                         %
% >> y = geospace(2, 7, 51);                                              %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Create an array with spacing increasing as a geometric series from 5 to %
% 32. Use 60 points.                                                      %
%                                                                         %
% >> y = geospace(5, 32, 60, 'reverse');                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [y] = geospace(x1, x2, n, varargin)


%% PARSE INPUTS

% Input defaults.
default.reverse = 0;
default.tol     = 1e-6;

% Input checks.
check.x1      = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric', 'sym'},                     ...
                     {'finite', 'real', 'scalar'});
check.x2      = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric', 'sym'},                     ...
                     {'finite', 'real', 'scalar'});
check.n       = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'integer', 'positive', 'scalar'});
check.reverse = @(x) any(validatestring(x, {'reverse'}));
check.tol     = @(x) validateattributes(x,                              ...
                     {'logical', 'numeric'},                            ...
                     {'finite', 'positive', 'real', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'x1'      ,                   check.x1      );
addRequired ( hParser , 'x2'      ,                   check.x2      );
addRequired ( hParser , 'n'       ,                   check.n       );
addOptional ( hParser , 'reverse' , default.reverse , check.reverse );
addParameter( hParser , 'tol'     , default.tol     , check.tol     );
parse(hParser, x1, x2, n, varargin{:});
clear check;

% Additional verifications.
nargoutchk(0,1);


%% CREATE THE GEOMETRIC SERIES ARRAY

% Define constants.
a0 = abs(x2 - x1);
a1 = 1 + a0;

% Strategy:
% * We need to solve the equation: r^n - a1*r + a0 = 0
% * I have thought of two, not necessarily guaranteed, methods. However,
%   they both seem to work quite well.
%   1) Use the roots() function of MATLAB and select the real root closest
%      to the common ratio for the infinite geometric series.
%      p = [1; zeros(n-2,1); -(1 + x2 - x1); (x2 - x1)];
%      r = roots(p);
%      r = min(abs(r(imag(r) == 0) - a0/a1));
%    * This method is not suitable for large array lengths due to the need
%      of finding all roots of a very large polynomial equation when we
%      really just need one.
%   2) Iteratively solve the equation by rearranging it as:
%      r_new = ((r_old)^n + a0)/a1
%    * This is the method I ultimately chose. I use the common ratio for
%      the infinite geometric series as an initial guess.

% Define the common ratio if the geometric series were infinite.
r0 = a0/a1;

% Determine the common ratio.
dr  = hParser.Results.tol + 1;
while dr > hParser.Results.tol
    r  = (r0^n + a0)/a1;
    dr = abs(r - r0);
    r0 = r;
end
clear a0 a1 dr r0;

% Determine the spacings (decreasing from x1 to x2).
dy = r.^((1:n-1));

% Evaluate the 'reverse' option.
if hParser.Results.reverse
    y = [0, cumsum(fliplr(dy))];
else
    y = [0, cumsum(dy)];
end

% Create the array and force the last point to be exact.
c    = (-1)^(double(x1) > double(x2));
y    = x1 + c*y.';
y(n) = x2;
clear c;


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
% 2022/03/07 -- (GDL) Added support for x1 > x2 and for vpa.              %
% 2022/03/03 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
