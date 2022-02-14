%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                                                                         %
% isgaussint                                                              %
% Data Interrogation                                                      %
% Examine data for genuine Gaussian integers                              %
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
% 2022/02/14 -- (GDL) Removed extra % found outside of preamble.          %
% 2021/06/06 -- (GDL) Added return false for empty arrays.                %
% 2021/06/04 -- (GDL) Added nargoutchk and future updates comments.       %
% 2021/05/27 -- (GDL) Added 'all' option for all() function.              %
% 2021/05/27 -- (GDL) Changed affiliation to ÉTS.                         %
% 2021/02/26 -- (GDL) Beta version of the code finalized.                 %
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
% val = isgaussint(S);                                                    %
% val = isgaussint(S, 'all');                                             %
% val = isgaussint(___, Name, Value);                                     %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Examine the elements of an input array and determine whether each is a  %
% genuine Gaussian integer (i.e., complex numbers whose real and          %
% imaginary parts are both integers). The 'all' option determines whether %
% the entire input array consists only of genuine Gaussian integers.      %
% Tolerances can be specified on the real and imaginary fractional parts  %
% of a number to account for typical small numerical errors. The default  %
% tolerances are zero. Elements of the input array that are either Inf or %
% NaN can be specified as genuine Gaussian integers or not. By default,   %
% such elements are considered to be genuine Gaussian integers.           %
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
%                for genuine Gaussian integers. Data types that are not   %
%                logical or numeric are accepted, but will return zero.   %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'all'          CASE-INSENSITIVE CHARACTER ARRAY                         %
%              ~ Specify 'all' to determine whether the entire input      %
%                array consists only of genuine Gaussian integers.        %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% 'inf'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether Inf values are to be treated as genuine  %
%                Gaussian integers (true, 1) or not (false, 0).           %
% ----------------------------------------------------------------------- %
% 'nan'          LOGICAL SCALAR                                           %
%                Default: true                                            %
%              ~ Specify whether NaN values are to be treated as genuine  %
%                Gaussian integers (true, 1) or not (false, 0).           %
% ----------------------------------------------------------------------- %
% 'tolIm'        NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance on the imaginary fractional part of a number.  %
% ----------------------------------------------------------------------- %
% 'tolRe'        NONNEGATIVE REAL SCALAR                                  %
%                Default: 0                                               %
%              ~ Tolerance on the real fractional part of a number.       %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'val'          LOGICAL N-DIMENSIONAL ARRAY/SCALAR                       %
%              ~ Output array. The elements of this array address whether %
%                the elements of the input array are genuine Gaussian     %
%                integers (1) or not (0). If the 'all' option is          %
%                specified, a logical scalar will instead address whether %
%                the entire input array consists only of genuine Gaussian %
%                integers (1) or not (0).                                 %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2i 0.71 NaN 10^-10+1i] %
% are genuine Gaussian integers. Use a tolerance of 10^-8 for both the    %
% real and imaginary fractional parts of a number and do not consider NaN %
% values as genuine Gaussian integers.                                    %
%                                                                         %
% >> x   = [2i 0.71 NaN 10^-10+1i];                                       %
% >> val = isgaussint(x, 'tolRe', 10^-8, 'tolIm', 10^-8, 'nan', 0);       %
% >> disp(val);                                                           %
%    1   0   0   1                                                        %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine whether the elements of the array x = [2+1i*10^-12 1+1i*20]   %
% are genuine Gaussian integers. Use a tolerance of 10^-8 for both the    %
% real and imaginary fractional parts of a number.                        %
%                                                                         %
% >> x   = [2+1i*10^-12 1+1i*20];                                         %
% >> val = isgaussint(x, 'tolRe', 10^-8, 'tolIm', 10^-8);                 %
% >> disp(val);                                                           %
%    1   1                                                                %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine whether the entire array x = [1i Inf+2i NaN*i 0] consists     %
% only of genuine Gaussian integers. Consider Inf and NaN values to be    %
% genuine Gaussian integers.                                              %
%                                                                         %
% >> x   = [1i Inf+2i NaN*i 0];                                           %
% >> val = isgaussint(x, 'all');                                          %
% >> disp(val);                                                           %
%    1                                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [val] = isgaussint(S, varargin)


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


%% GAUSSIAN INTEGER INTERROGATION

% If the input is not logical nor numeric, or if the input is empty, then
% it cannot be a Gaussian integer. The char data type is not accepted as an
% integer by MATfluids.
if ~(islogical(S) || isnumeric(S)) || isempty(S)
    val = 0;
    return;
end

% Initialize a cell to hold the conditions.
C = cell(6,1);

% Re: (Allow Inf?) AND (Is Inf?)
C{1} =  hParser.Results.inf                                             ...
     &  isinf(real(S));
% Re: (Allow NaN?) AND (Is NaN?)
C{2} =  hParser.Results.nan                                             ...
     &  isnan(real(S));
% Re: (Not Inf?) AND (Not NaN?) AND (No Mod?)
C{3} = ~isinf(real(S))                                                  ...
     & ~isnan(real(S))                                                  ...
     & (abs(mod(real(S),1)) <= hParser.Results.tolRe);
% Im: (Allow Inf?) AND (Is Inf?)
C{4} =  hParser.Results.inf                                             ...
     &  isinf(imag(S));
% Im: (Allow NaN?) AND (Is NaN?)
C{5} =  hParser.Results.nan                                             ...
     &  isnan(imag(S));
% Im: (Not Inf?) AND (Not NaN?) AND (No Mod?)
C{6} = ~isinf(imag(S))                                                  ...
     & ~isnan(imag(S))                                                  ...
     & (abs(mod(imag(S),1)) <= hParser.Results.tolIm);

% (Inf?) OR (NaN?) OR (Integer?), Re AND Im
val = (C{1} | C{2} | C{3}) & (C{4} | C{5} | C{6});
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
