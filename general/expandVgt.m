%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% expandVgt                                                               %
% MATfluids Structure Interrogation                                       %
% Expand contents of velocity gradient structure to most general form     %
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
% evgt = expandVgt(vgt);                                                  %
% evgt = expandVgt(vgt, val);                                             %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Expand a valid MATfluids velocity gradient structure to include all     %
% possible fields ["ux" "uy" "uz" "vx" "vy" "vz" "wx" "wy" "wz"]. Missing %
% fields are added as empty arrays by default, though the user can also   %
% specify a scalar value.                                                 %
%                                                                         %
% Compatibility:                                                          %
% MATLAB R2019b or later.                                                 %
%                                                                         %
% Dependencies:                                                           %
% examineVgt                                                              %
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
% 'vgt'          STRUCT ARRAY (1 X 1)                                     %
%              ~ Velocity gradient tensor. The structure array may        %
%                contain the following fields:                            %
%                1) ux : Gradient of u component in the x direction       %
%                2) uy : Gradient of u component in the y direction       %
%                3) uz : Gradient of u component in the z direction       %
%                4) vx : Gradient of v component in the x direction       %
%                5) vy : Gradient of v component in the y direction       %
%                6) vz : Gradient of v component in the z direction       %
%                7) wx : Gradient of w component in the x direction       %
%                8) wy : Gradient of w component in the y direction       %
%                9) wz : Gradient of w component in the z direction       %
%                Each field is a one to four dimensional array, in ndgrid %
%                format (i.e., dimensions represent [x y z t]). All       %
%                constituent arrays must have the same size. Not all      %
%                fields need be present in the velocity gradient          %
%                structure, however the order must follow [ux uy uz vx vy %
%                vz wx wy wz]. For example, the velocity gradient         %
%                structure may contain fields ordered as [ux uy vx vy]    %
%                but not as [vx vy ux uy]. Empty, scalar and complex      %
%                arrays are permitted.                                    %
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
% 'evgt'         STRUCT ARRAY (1 X 1)                                     %
%              ~ Expanded velocity gradient tensor. The structure array   %
%                contains all the following fields:                       %
%                1) ux : Gradient of u component in the x direction       %
%                2) uy : Gradient of u component in the y direction       %
%                3) uz : Gradient of u component in the z direction       %
%                4) vx : Gradient of v component in the x direction       %
%                5) vy : Gradient of v component in the y direction       %
%                6) vz : Gradient of v component in the z direction       %
%                7) wx : Gradient of w component in the x direction       %
%                8) wy : Gradient of w component in the y direction       %
%                9) wz : Gradient of w component in the z direction       %
%                Fields available in the input velocity gradient          %
%                structure are copied whereas fields that aren't are      %
%                added as empty arrays.                                   %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Expand a MATfluids velocity gradient structure initially containing     %
% only the ["ux" "uy" "vx" "vy"] fields.                                  %
%                                                                         %
% >> vgt.ux = 1;                                                          %
% >> vgt.uy = [1 2 3; 4 5 6; 7 8 9];                                      %
% >> vgt.vx = [1 2 1; 2 1 2; 1 2 1];                                      %
% >> vgt.vy = [];                                                         %
% >> vgt = expandVgt(vgt);                                                %
% >> disp(vgt);                                                           %
%    ux: 1                                                                %
%    uy: [3×3 double]                                                     %
%    uz: []                                                               %
%    vx: [3×3 double]                                                     %
%    vy: []                                                               %
%    vz: []                                                               %
%    wx: []                                                               %
%    wy: []                                                               %
%    wz: []                                                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [evgt] = expandVgt(vgt, varargin)


%% PARSE INPUTS

% Input defaults.
default.val = [];

% Input checks.
check.vgt = @(x) examineVgt(x);
check.val = @(x) validateattributes(x,                                  ...
                 {'logical', 'numeric'},                                ...
                 {'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser , 'vgt' ,              check.vgt );
addOptional ( hParser , 'val' , default.val, check.val );
parse(hParser, vgt, varargin{:});
clear check;

% Additional verifications.
narginchk(1,2);
nargoutchk(0,1);


%% INITIALIZATIONS

evgt = vgt;                 % Expanded velocity gradient structure
val  = hParser.Results.val; % Assigned value to missing fields


%% EXPAND VELOCITY GRADIENT STRUCTURE

% Add missing fields.
if ~isfield(evgt, 'ux'), evgt.ux = val; end
if ~isfield(evgt, 'uy'), evgt.uy = val; end
if ~isfield(evgt, 'uz'), evgt.uz = val; end
if ~isfield(evgt, 'vx'), evgt.vx = val; end
if ~isfield(evgt, 'vy'), evgt.vy = val; end
if ~isfield(evgt, 'vz'), evgt.vz = val; end
if ~isfield(evgt, 'wx'), evgt.wx = val; end
if ~isfield(evgt, 'wy'), evgt.wy = val; end
if ~isfield(evgt, 'wz'), evgt.wz = val; end
clear val;

% Order the fields in MATfluids format.
evgt = orderfields(evgt, {'ux', 'uy', 'uz', 'vx', 'vy', 'vz', 'wx',     ...
                          'wy', 'wz'});


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
