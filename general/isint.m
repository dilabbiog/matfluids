%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isint                                                                   %
% Data Interrogation                                                      %
% Examine data for genuine integers                                       %
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
% val = isint(S);                                                         %
% val = isint(S, 'all');                                                  %
% val = isint(___, Name, Value);                                          %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine the elements of an input array and determine whether each is a  %
% genuine integer. The 'all' option determines whether the entire input   %
% array consists only of genuine integers. A tolerance can be specified   %
% on the fractional part of a real number to account for typical small    %
% numerical errors. Likewise, a tolerance can be specified on the         %
% imaginary part of a number. The default tolerances are zero. Elements   %
% of the input array that are either Inf or NaN can be specified as       %
% genuine integers or not. By default, such elements are considered to be %
% genuine integers.                                                       %
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
% 'S'            LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. The elements of this array will be examined %
%                for genuine integers. Data types that are not logical or %
%                numeric are accepted, but will return zero.              %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether the entire input      %
%                array consists only of genuine integers.                 %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'inf'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether Inf values are to be treated as genuine  %
%                integers (true, 1) or not (false, 0).                    %
% ----------------------------------------------------------------------- %
% 'nan'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether NaN values are to be treated as genuine  %
%                integers (true, 1) or not (false, 0).                    %
% ----------------------------------------------------------------------- %
% 'tolIm'        NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance on the imaginary part of a complex number.     %
% ----------------------------------------------------------------------- %
% 'tolRe'        NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance on the fractional part of a real number.       %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL N-DIMENSIONAL ARRAY/SCALAR                       %
%              ~ Output array. The elements of this array address whether %
%                the elements of the input array are genuine integers (1) %
%                or not (0). If the 'all' option is specified, a logical  %
%                scalar will instead address whether the entire input     %
%                array consists only of genuine integers (1) or not (0).  %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2 Inf 0.1 2+10^-12]    %
% are genuine integers. Use a tolerance of 10^-8 for the fractional parts %
% of real numbers and do not consider Inf values as genuine integers.     %
%                                                                         %
% >> x   = [2 Inf 0.1 2+10^-12];                                          %
% >> val = isint(x, 'tolRe', 10^-8, 'inf', 0);                            %
% >> disp(val);                                                           %
%    1   0   0   1                                                        %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2+1i*10^-12 1+1i*20]   %
% are genuine integers. Use a tolerance of 10^-8 for the imaginary parts  %
% of complex numbers.                                                     %
%                                                                         %
% >> x   = [2+1i*10^-12 1+1i*20];                                         %
% >> val = isint(x, 'tolIm', 10^-8);                                      %
% >> disp(val);                                                           %
%    1   0                                                                %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the entire array x = [1 -2 NaN 0] consists only of    %
% genuine integers. Consider NaN values to be genuine integers.           %
%                                                                         %
% >> x   = [1 -2 NaN 0];                                                  %
% >> val = isint(x, 'all');                                               %
% >> disp(val);                                                           %
%    1                                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = isint(S, varargin)


%% PARSE INPUTS

% Input defaults.
default.all   = 0;
default.inf   = 1;
default.nan   = 1;
default.tolIm = 0;
default.tolRe = 0;

% Input checks.
check.all   = @(x) any(validatestring(x, {'all'}));
check.inf   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'binary', 'scalar'});
check.nan   = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'binary', 'scalar'});
check.tolIm = @(x) validateattributes(x,                                ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'nonnegative', 'real', 'scalar'}) ;
check.tolRe = @(x) validateattributes(real(x),                          ...
                   {'logical', 'numeric'},                              ...
                   {'finite', 'nonnegative', 'real', 'scalar'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'S'                                   );
addOptional ( hParser, 'all'   , default.all   , check.all   );
addParameter( hParser, 'inf'   , default.inf   , check.inf   );
addParameter( hParser, 'nan'   , default.nan   , check.nan   );
addParameter( hParser, 'tolIm' , default.tolIm , check.tolIm );
addParameter( hParser, 'tolRe' , default.tolRe , check.tolRe );
parse(hParser, S, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,1);


%% INTEGER INTERROGATION

% If the input is not logical nor numeric, or if the input is empty, then
% it cannot be an integer. The char data type is not accepted as an integer
% by MATfluids.
if ~(islogical(S) || isnumeric(S)) || isempty(S)
    val = 0;
    return;
end

% Initialize a cell to hold the conditions.
C = cell(3,1);

% (Allow Inf?) AND (Is Inf?) AND (No Im?)
C{1} =  hParser.Results.inf                                             ...
     &  isinf(real(S))                                                  ...
     & ~isinf(imag(S))                                                  ...
     & (abs(imag(S)) <= hParser.Results.tolIm);
% (Allow NaN?) AND (Is NaN?) AND (No Im?)
C{2} =  hParser.Results.nan                                             ...
     &  isnan(real(S))                                                  ...
     & ~isnan(imag(S))                                                  ...
     & (abs(imag(S)) <= hParser.Results.tolIm);
% (Not Inf?) AND (Not NaN?) AND (No Mod?) AND (No Im?)
C{3} = ~isinf(real(S))                                                  ...
     & ~isnan(real(S))                                                  ...
     & (abs(mod(real(S),1)) <= hParser.Results.tolRe)                   ...
     & (abs(imag(S)) <= hParser.Results.tolIm);

% (Inf?) OR (NaN?) OR (Integer?)
val = C{1} | C{2} | C{3};
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
% 2021/06/06 -- (GDL) Added return false for empty arrays.                %
% 2021/06/04 -- (GDL) Added nargoutchk and future updates comments.       %
% 2021/05/27 -- (GDL) Added 'all' option for all() function.              %
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
