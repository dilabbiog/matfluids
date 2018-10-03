%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%          EXTRACT A SPECIFIED FIELD FROM THE WORKSPACE VARIABLES         %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 3rd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2018 Giuseppe Di Labbio                                   %
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
% WSprop = extractWSfield(property);                                      %
% WSprop = extractWSfield(property, workspace);                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function extracts a specific property ('name', 'size', 'bytes',    %
% 'class', 'global', 'sparse', 'complex', 'nesting', 'persistent') from   %
% the variables in the MATLAB workspace. The user has the added option of %
% using MATLAB's base workspace or the workspace of the caller function   %
% or script. The function is equivalent to the commands                   %
% >> tmp = whos;                                                          %
% >> tmp = extractfield(tmp, property).';                                 %
% where extractfield() is a function in MATLAB's Mapping Toolbox [1].     %
% Here, there is no dependence on MATLAB's Mapping Toolbox whatsoever.    %
%                                                                         %
% References:                                                             %
% [1] The MathWorks, Inc. (2016). Functions - Alphabetical List. In       %
%     Mapping Toolbox Reference: MATLAB R2016a (297-298). Natick, MA: The %
%     MathWorks, Inc.                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'property'   - STRING                                                   %
%              - Field name of the property to extract from the MATLAB    %
%                workspace in question. The user must select from:        %
%               'bytes', 'class', 'complex', 'global', 'name', 'nesting', %
%               'persistent', 'size', 'sparse'                            %
% ----------------------------------------------------------------------- %
% 'workspace'  - STRING                                                   %
%              - Workspace frome which to extract the specified property. %
%                The user can select either 'base' or 'caller' to use     %
%                MATLAB's base workspace or the caller function/script's  %
%                workspace respectively.                                  %
%              - Default: 'caller'                                        %
% ----------------------------------------------------------------------- %
% 'WSprop'     - 1D CELL ARRAY                                            %
%              - One-dimensional cell array holding the extracted         %
%                property of the variables from the selected workspace.   %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Create a few variables and extract the names and memory size of all the %
% variables in the workspace, including the variables that will be        %
% created to hold the names and memory size. Additionally, express the    %
% size of the variables in kB.                                            %
%                                                                         %
% >> a = ones(5,2);                                                       %
% >> x = 'Hello World!';                                                  %
% >> f = @(x) 3*x + 2;                                                    %
% >> q = a - f(12);                                                       %
% >> names = sort([extractWSfield('name'); 'names'; 'mem']);              %
% >> mem = 0;                                                             %
% >> mem = cell2mat(extractWSfield('bytes'))/1000;                        %
% >> mem = cell2mat(extractWSfield('bytes'))/1000;                        %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extractWSfield
function [var] = extractWSfield(fieldname, varargin)

%% Parse the inputs.
if ~nargin
    error('InputError:noInput', 'Not enough input arguments.');
elseif nargin == 1
    if ischar(fieldname)
        chk = { 'name' 'size' 'bytes' 'class' 'global' 'sparse'         ...
                'complex' 'nesting' 'persistent'                }.';
        if ~max(strcmpi(fieldname, chk))
            error('InputError:notField',                                ...
                 ['The field name could only be one of the following: ' ...
                  'name, size, bytes, class, global, sparse, complex, ' ...
                  'nesting or persistent.']);
        end
    else
        error('InputError:notString1', 'The field name must be a string.');
    end
elseif nargin == 2
    if ischar(fieldname)
        chk = { 'name' 'size' 'bytes' 'class' 'global' 'sparse'         ...
                'complex' 'nesting' 'persistent'                }.';
        if ~max(strcmpi(fieldname, chk))
            error('InputError:notField',                                ...
                 ['The field name could only be one of the following: ' ...
                  'name, size, bytes, class, global, sparse, complex, ' ...
                  'nesting or persistent.']);
        end
    else
        error('InputError:notString1', 'The field name must be a string.');
    end
    if ischar(varargin{1})
        chk = {'base' 'caller'}.';
        if ~max(strcmpi(varargin{1}, chk))
            error('InputError:notWorkspace',                            ...
                 ['There are only two workspaces to choose from: base ' ...
                  'or caller.']);
        end
    else
        error('InputError:notString2', 'The workspace must be a string.');
    end
else
    error('InputError:tooMany', 'Too many input arguments.');
end

%% Determine which workspace to extract the selected field from.
useBase   = 0;
useCaller = 0;
if nargin == 2
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'caller')
        useCaller = 1;
    elseif ischar(varargin{1}) && strcmpi(varargin{1}, 'base')
        useBase = 1;
    end
end

%% Extract the selected field.
if nargin == 1 || (nargin == 2 && useCaller)
    tmp = evalin('caller', 'whos');
    var = cell(length(tmp),1);
    for k = 1:length(tmp)
        var{k} = tmp(k).(fieldname);
    end
elseif nargin == 2 && useBase
    tmp = evalin('base', 'whos');
    var = cell(length(tmp),1);
    for k = 1:length(tmp)
        var{k} = tmp(k).(fieldname);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A.
% Reason(s)
% * N/A.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) 164, 170                                                        %
% * Following MATLAB's recommendation: "Use dynamic fieldnames with       %
%   structures instead of GETFIELD", lines 164 and 170 were selected over %
%   using getfield(tmp(k), fieldname).                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
