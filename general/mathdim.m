%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% mathdim                                                                 %
% Array Properties                                                        %
% Determine the mathematical dimension of an array                        %
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
% CHANGE LOG                                                              %
%                                                                         %
% 2021/05/27 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2021 Giuseppe Di Labbio                                   %
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
% dim = mathdim(A);                                                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Determine the true mathematical dimension of an array. In other words,  %
% a scalar will have a dimension of 0, a vector will have a dimension of  %
% 1 (regardless of which direction its length lies), a matrix will have a %
% dimension of 2 and so on. For arrays of dimension 2 or greater, the     %
% mathematical dimension is taken as ndims(squeeze()).                    %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R20019b or later.                                                %
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
% 'A'            LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. The mathematical dimension of this array    %
%                will be determined.                                      %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'dim'          NONNEGATIVE INTEGER SCALAR                               %
%              ~ Output scalar. True mathematical dimension of a logical  %
%                or numeric array, being 0 for a scalar, 1 for a vector,  %
%                2 for a matrix and so on.                                %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine the true mathematical dimension of the scalar x = 20.         %
%                                                                         %
% >> x   = 20;                                                            %
% >> dim = mathdim(x);                                                    %
% >> disp(dim);                                                           %
%      0                                                                  %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine the true mathematical dimension of the array x = [1 2 3].     %
%                                                                         %
% >> x   = [1 2 3];                                                       %
% >> dim = mathdim(x);                                                    %
% >> disp(dim);                                                           %
%      1                                                                  %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine the true mathematical dimension of the three-dimensional      %
% array x = permute([1 2 3], [3 1 2]).                                    %
%                                                                         %
% >> x   = permute([1 2 3], [3 1 2]);                                     %
% >> dim = mathdim(x);                                                    %
% >> disp(dim);                                                           %
%      1                                                                  %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 4                                                               %
%                                                                         %
% Determine the true mathematical dimension of the four-dimensional array %
% x = permute([1 2 3; 4 5 6; 7 8 9], [4 3 1 2]).                          %
%                                                                         %
% >> x   = permute([1 2 3; 4 5 6; 7 8 9], [4 3 1 2]);                     %
% >> dim = mathdim(x);                                                    %
% >> disp(dim);                                                           %
%      2                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dim] = mathdim(A)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
check.A = @(x) validateattributes(x,                                    ...
               {'logical', 'numeric'},                                  ...
               {});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'A', check.A );
parse(hParser, A);
clear check;

% Additional verifications.
% N/A


%% DETERMINE MATHEMATICAL DIMENSION

if     isscalar(A),          dim = 0;
elseif isvector(squeeze(A)), dim = 1;
else,                        dim = ndims(squeeze(A));
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
