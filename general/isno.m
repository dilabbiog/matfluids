%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isno                                                                    %
% Data Interrogation                                                      %
% Examine data for a form of "no"                                         %
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
% val = isno(S);                                                          %
% val = isno(S, 'all');                                                   %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine whether a character or string input is a form of "no". The      %
% 'all' option determines whether an entire input string array consists   %
% only of forms of "no". In the case of a character array, only simple    %
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
%                for a form of "no". The following words will return      %
%                true:                                                    %
%                'Laa', 'N', 'Na', 'Nah', 'Nao', 'Nay', 'Ne', 'Nee',      %
%                'Neen', 'Negative', 'Nei', 'Nein', 'Nej', 'Nem', 'Nie',  %
%                'No', 'Non', Nope'                                       %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether the entire input      %
%                string array consists only of forms of "no".             %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL SCALAR                                           %
%              - Logical answer to whether the interrogated character     %
%                array or string is a form of "no" (1) or not (0).        %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Add two numbers if a governing variable is set to a form of "no".       %
%                                                                         %
% >> someCond = 'N';                                                      %
% >> if isno(someCond)                                                    %
%        disp(2 + 5);                                                     %
%    end                                                                  %
%      7                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = isno(S, varargin)


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
verify.char = ischar(S) & isrow(S);


%% YES INTERROGATION

% If the input is not character (row vector) or string, then it cannot be a
% form of "no".
if ~(verify.char || isstring(S))
    val = 0;
    return;
end
clear verify;

% Define words that qualify as a form of "no".
C = {'Laa', 'N', 'Na', 'Nah', 'Nao', 'Nay', 'Ne', 'Nee', 'Neen',        ...
     'Negative', 'Nei', 'Nein', 'Nej', 'Nem', 'Nie', 'No', 'Non',       ...
     'Nope'}.';

% Examine whether the input qualifies as a form of "no".
val = ismember(lower(S), lower(C));
clear C;

% Evaluate the 'all' option.
if hParser.Results.all
    val = all(val, 'all');
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
