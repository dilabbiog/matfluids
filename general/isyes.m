%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isyes                                                                   %
% Data Interrogation                                                      %
% Examine data for a form of "yes"                                        %
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
% val = isyes(S);                                                         %
% val = isyes(S, 'all');                                                  %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a character or string input is a form of "yes". The     %
% 'all' option determines whether an entire input string array consists   %
% only of forms of "yes". In the case of a character array, only simple   %
% words (row vector character arrays) are considered.                     %
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
% 'S'            CASE-INSENSITIVE CHARACTER/STRING ARRAY                  %
%              ~ Input array. The elements of this array will be examined %
%                for a form of "yes". The following words will return     %
%                true:                                                    %
%                'Baleh', 'Da', 'Etiam', 'Go', 'Hai', 'Ja', 'K', 'Naam',  %
%                'Ok', 'Okay', 'Oui', 'Please', 'Positive', 'Si', 'Sim',  %
%                'Sure', 'True', 'Y', 'Ya', 'Yay', 'Ye', 'Yeah', 'Yeh',   %
%                'Yes', 'Yup'                                             %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether the entire input      %
%                string array consists only of forms of "yes".            %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL SCALAR                                           %
%              ~ Logical answer to whether the interrogated character     %
%                array or string is a form of "yes" (1) or not (0).       %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Add two numbers if a governing variable is set to a form of "yes".      %
%                                                                         %
% >> someCond = 'Y';                                                      %
% >> if isyes(someCond)                                                   %
%        disp(2 + 5);                                                     %
%    end                                                                  %
%      7                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = isyes(S, varargin)


%% PARSE INPUTS

% Input defaults.
default.all = 0;

% Input checks.
check.all = @(x) any(validatestring(x, {'all'}));

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'S'                             );
addOptional ( hParser, 'all' , default.all , check.all );
parse(hParser, S, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);
verify.char = ischar(S) & isrow(S);


%% YES INTERROGATION

% If the input is not character (row vector) or string, then it cannot be a
% form of "yes".
if ~(verify.char || isstring(S))
    val = 0;
    return;
end
clear verify;

% Define words that qualify as a form of "yes".
C = {'Baleh', 'Da', 'Etiam', 'Go', 'Hai', 'Ja', 'K', 'Naam', 'Ok',      ...
     'Okay', 'Oui', 'Please', 'Positive', 'Si', 'Sim', 'Sure', 'True',  ...
     'True', 'Y', 'Ya', 'Yay', 'Ye', 'Yeah', 'Yeh', 'Yes', 'Yup'}.';

% Examine whether the input qualifies as a form of "yes".
val = ismember(lower(S), lower(C));
clear C;

% Evaluate the 'all' option.
if hParser.Results.all
    val = all(val, 'all');
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
%#ok<*N/A>                                                                %
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
% 2022/02/22 -- (GDL) Moved change log and future updates to bottom,      %
%                     reformatted notes.                                  %
% 2021/06/04 -- (GDL) Added nargoutchk and future updates comments.       %
% 2021/05/27 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
