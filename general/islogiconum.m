%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% islogiconum                                                             %
% Data Type Interrogation                                                 %
% Examine whether data type is logical or numeric                         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% D�partement de g�nie m�canique                                          %
% �cole de technologie sup�rieure (�TS)                                   %
% Montr�al, Qu�bec                                                        %
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
% val = islogiconum(S);                                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether an input has the logical or numeric data type. Returns  %
% true if either islogical() or isnumeric() return true.                  %
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
% 'S'            ANY INPUT                                                %
%              ~ Input variable. The data type will be examined for       %
%                whether it is logical or numeric.                        %
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
% 'val'          LOGICAL SCALAR                                           %
%              ~ Output logical scalar. The value will be 1 (true) if the %
%                input data type is either logical or numeric and 0       %
%                (false) otherwise.                                       %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the array x = [2 Inf NaN 10^-12 -Inf 2*i -3.1] is     %
% logical or numeric.                                                     %
%                                                                         %
% >> x   = [2 Inf NaN 10^-12 -Inf 2*i -3.1];                              %
% >> val = islogiconum(x);                                                %
% >> disp(val);                                                           %
%    1                                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the logical array x = logical([0 true 1 false]) is    %
% logical or numeric.                                                     %
%                                                                         %
% >> x   = logical([0 true 1 false]);                                     %
% >> val = islogiconum(x);                                                %
% >> disp(val);                                                           %
%    1                                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the cell array x = {1 NaN 2.4} is logical or numeric. %
%                                                                         %
% >> x   = {1 NaN 2.4};                                                   %
% >> val = islogiconum(x);                                                %
% >> disp(val);                                                           %
%    0                                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = islogiconum(S)


%% PARSE INPUTS

% Input defaults.
% N/A

% Input checks.
% N/A

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'S' );
parse(hParser, S);

% Additional verifications.
nargoutchk(0,1);


%% LOGICONUMERIC INTERROGATION

% (Is numeric?) OR (Is logical?)
val = isnumeric(S) | islogical(S);


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
% 2022/02/23 -- (GDL) Removed message suppression in file, prefer line.   %
% 2022/02/22 -- (GDL) Moved change log and future updates to bottom,      %
%                     reformatted notes.                                  %
% 2021/06/04 -- (GDL) Added nargoutchk and future updates comments.       %
% 2021/05/27 -- (GDL) Added simple input parser.                          %
% 2021/05/27 -- (GDL) Changed affiliation to �TS.                         %
% 2021/02/26 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
