%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isnegreal                                                               %
% Data Interrogation                                                      %
% Examine data for genuine negative real numbers                          %
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
% CHANGE LOG                                                              %
%                                                                         %
% 2021/06/06 -- (GDL) Added return false for empty arrays.                %
% 2021/06/04 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
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
% val = isnegreal(S);                                                     %
% val = isnegreal(S, 'all');                                              %
% val = isnegreal(___, Name, Value);                                      %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine the elements of an input array and determine whether each is a  %
% genuine negative real number. The 'all' option determines whether the   %
% entire input array consists only of genuine negative real numbers. A    %
% tolerance can be specified on the fractional part of a real number to   %
% account for typical small numerical errors (i.e., what qualifies as     %
% zero). Likewise, a tolerance can be specified on the imaginary part of  %
% a number. The default tolerances are zero. Elements of the input array  %
% that are either -Inf or NaN can be specified as genuine negative real   %
% numbers or not. By default, such elements are considered to be genuine  %
% negative real numbers.                                                  %
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
%                for genuine negative real numbers. Data types that are   %
%                not logical or numeric are accepted, but will return     %
%                zero.                                                    %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether the entire input      %
%                array consists only of genuine negative real numbers.    %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'inf'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether -Inf values are to be treated as genuine %
%                negative real numbers (true, 1) or not (false, 0).       %
% ----------------------------------------------------------------------- %
% 'nan'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether NaN values are to be treated as genuine  %
%                negative real numbers (true, 1) or not (false, 0).       %
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
%                the elements of the input array are genuine negative     %
%                real numbers (1) or not (0). If the 'all' option is      %
%                specified, a logical scalar will instead address whether %
%                the entire input array consists only of genuine negative %
%                real numbers (1) or not (0).                             %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2 Inf 0.1 1e-12 -Inf]  %
% are genuine negative real numbers. Consider numbers below 1e-8 as zero. %
%                                                                         %
% >> x   = [2 Inf 0.1 1e-12 -Inf];                                        %
% >> val = isnegreal(x, 'tolRe', 1e-8);                                   %
% >> disp(val);                                                           %
%    0   0   0   0   1                                                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2+1i*1e-12 -Inf 1+1i]  %
% are genuine negative real numbers. Use a tolerance of 1e-8 for the      %
% imaginary parts of complex numbers and do not consider -Inf values to   %
% be genuine negative real numbers.                                       %
%                                                                         %
% >> x   = [2+1i*1e-12 -Inf 1+1i];                                        %
% >> val = isnegreal(x, 'tolIm', 1e-8, 'inf', 0);                         %
% >> disp(val);                                                           %
%    0   0   0                                                            %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the entire array x = [-1.3 -2 NaN -0.7 0] contains    %
% onlybgenuine negative real numbers. Consider NaN values to be genuine   %
% negative real numbers.                                                  %
%                                                                         %
% >> x   = [-1.3 -2 NaN -0.7 0];                                          %
% >> val = isnegreal(x, 'all');                                           %
% >> disp(val);                                                           %
%    0                                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = isnegreal(S, varargin)


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


%% NEGATIVE REAL NUMBER INTERROGATION

% If the input is not logical nor numeric, or if the input is empty, then
% it cannot be a negative real number. The char data type is not accepted
% as an integer by MATfluids.
if ~(islogical(S) || isnumeric(S)) || isempty(S)
    val = 0;
    return;
end

% Initialize a cell to hold the conditions.
C = cell(3,1);

% (Allow Inf?) AND (Is Inf?) AND (Is -?) AND (No Im?)
C{1} =  hParser.Results.inf                                             ...
     &  isinf(real(S))                                                  ...
     & (sign(real(S)) == -1)                                            ...
     & ~isinf(imag(S))                                                  ...
     & (abs(imag(S)) <= hParser.Results.tolIm);
% (Allow NaN?) AND (Is NaN?) AND (No Im?)
C{2} =  hParser.Results.nan                                             ...
     &  isnan(real(S))                                                  ...
     & ~isnan(imag(S))                                                  ...
     & (abs(imag(S)) <= hParser.Results.tolIm);
% (Not Inf?) AND (Not NaN?) AND (Is -?) AND (No Im?)
C{3} = ~isinf(real(S))                                                  ...
     & ~isnan(real(S))                                                  ...
     & (real(S) < -hParser.Results.tolRe)                               ...
     & (abs(imag(S)) <= hParser.Results.tolIm);

% (-Inf?) OR (NaN?) OR (Negative Real?)
val = C{1} | C{2} | C{3};
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
% * N/A                                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%