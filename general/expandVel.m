%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% expandVel                                                               %
% MATfluids Structure Interrogation                                       %
% Expand contents of velocity structure to most general form              %
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
% Copyright (C) 2023 Giuseppe Di Labbio                                   %
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
% evel = expandVel(vel);                                                  %
% evel = expandVel(vel, val);                                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Expand a valid MATfluids velocity structure to include all possible     %
% fields ["u" "v" "w"]. Missing fields are added as empty arrays by       %
% default, though the user can also specify a scalar value.               %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% examineVel                                                              %
% islogiconum                                                             %
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
% 'vel'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity field. The structure array may contain the      %
%                following fields:                                        %
%                1) u : Velocity component in the x direction             %
%                2) v : Velocity component in the y direction             %
%                3) w : Velocity component in the z direction             %
%                Each field is a one to four dimensional array, in ndgrid %
%                format (i.e., dimensions represent [x y z t]). All       %
%                constituent arrays must have the same size. Not all      %
%                fields need be present in the velocity structure,        %
%                however the order must follow [u v w]. For example, the  %
%                velocity structure may contain fields ordered as [u w]   %
%                but not as [w u]. Empty, scalar and complex arrays are   %
%                permitted.                                               %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'val'          SCALAR                                                   %
%              ~ Input scalar. Missing fields of the velocity structure   %
%                will be assigned this value.                             %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'evel'         STRUCT ARRAY (1 X 1)                                     %
%              ~ Expanded velocity field. The structure array contains    %
%                all the following fields:                                %
%                1) u : Velocity component in the x direction             %
%                2) v : Velocity component in the y direction             %
%                3) w : Velocity component in the z direction             %
%                Fields available in the input velocity structure are     %
%                copied whereas fields that aren't are added as empty     %
%                arrays.                                                  %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Expand a MATfluids velocity structure initially containing only the     %
% ["u" "v"] fields.                                                       %
%                                                                         %
% >> vel.u = 1;                                                           %
% >> vel.v = [1 2 3; 4 5 6; 7 8 9];                                       %
% >> vel = expandVel(vel);                                                %
% >> disp(vel);                                                           %
%    u: 1                                                                 %
%    v: [3x3 double]                                                      %
%    w: []                                                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [evel] = expandVel(vel, varargin)


%% PARSE INPUTS

% Input defaults.
default.val = [];

% Input checks.
check.vel = @(x) examineVel(x);
check.val = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'vel' ,              check.vel );
addOptional ( hParser , 'val' , default.val, check.val );
parse(hParser, vel, varargin{:});
clear check;

% Additional verifications.
narginchk(1,2);
nargoutchk(0,1);


%% INITIALIZATIONS

evel = vel;                 % Expanded velocity structure
val  = hParser.Results.val; % Assigned value to missing fields


%% EXPAND VELOCITY STRUCTURE

% Add missing fields.
if ~isfield(evel, 'u'), evel.u = val; end
if ~isfield(evel, 'v'), evel.v = val; end
if ~isfield(evel, 'w'), evel.w = val; end
clear val;

% Order the fields in MATfluids format.
evel = orderfields(evel, {'u', 'v', 'w'});


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
% 2023/09/22 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
