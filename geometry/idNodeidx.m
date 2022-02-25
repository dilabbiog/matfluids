%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            GEOMETRY TOOLBOX                             %
%                                                                         %
% idNodeidx                                                               %
% Node Property Identification                                            %
% Examine geometry for indices of leading, trailing and interior nodes    %
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
% [idxl, idxi, idxt] = idNodeidx(geom);                                   %
% [idxl, idxi, idxt] = idNodeidx(geom, dir);                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Determine the linear indices of the leading, interior and trailing      %
% nodes of an input array. By default, the indices are listed in the      %
% order they appear along the first dimension (since this is how linear   %
% indices are counted by MATLAB). Alternatively, the indices may be       %
% listed in the order they appear along a different direction.            %
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
% 'geom'         LOGICAL/NUMERIC N-DIMENSIONAL ARRAY                      %
%              ~ Input array. Geometry or mask. Values are 1 within the   %
%                geometry and 0 outside the geometry.                     %
% ======================================================================= %
% Input Arguments (Optional):                                             %
% ----------------------------------------------------------------------- %
% 'dir'          POSITIVE INTEGER SCALAR                                  %
%                Default: 1                                               %
%              ~ Input scalar. Direction along which to list the indices. %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'idxl'         NUMERIC COLUMN VECTOR                                    %
%              ~ Output array. List of linear indices of leading nodes in %
%                the order they appear along a specified direction.       %
% ----------------------------------------------------------------------- %
% 'idxi'         NUMERIC COLUMN VECTOR                                    %
%              ~ Output array. List of linear indices of interior nodes   %
%                in the order they appear along a specified direction.    %
% ----------------------------------------------------------------------- %
% 'idxt'         NUMERIC COLUMN VECTOR                                    %
%              ~ Output array. List of linear indices of trailing nodes   %
%                in the order they appear along a specified direction.    %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in one dimension.                                               %
%                                                                         %
% >> geom               = [ones(20,1); zeros(10,1); ones(6,1)];           %
% >> [idxl, idxi, idxt] = idNodeidx(geom);                                %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in two dimensions.                                              %
%                                                                         %
% >> geom               = zeros(51, 101);                                 %
% >> geom(5:20,5:20)    = 1;                                              %
% >> geom(40:45,55:80)  = 1;                                              %
% >> [idxl, idxi, idxt] = idNodeidx(geom, 1);                             %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in two dimensions. List in the indices in the order they appear %
% along the second dimension.                                             %
%                                                                         %
% >> geom               = zeros(51, 101);                                 %
% >> geom(5:20,5:20)    = 1;                                              %
% >> geom(40:45,55:80)  = 1;                                              %
% >> [idxl, idxi, idxt] = idNodeidx(geom, 2);                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [idxl, idxi, idxt] = idNodeidx(geom, varargin)


%% PARSE INPUTS

% Input defaults.
default.dir = 1;

% Input checks.
check.dir  = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'positive', 'scalar'});
check.geom = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty', 'binary'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'geom' , check.geom              );
addOptional ( hParser, 'dir'  , default.dir , check.dir );
parse(hParser, geom, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,3);


%% INDEX SEARCH

% Determine the size of the array.
sz = size(geom);

% Set the search direction.
dir = hParser.Results.dir;

% Permute the geometry array.
pm = 1:ndims(geom);
if dir > 1
    pm(dir) = [];
    pm      = [dir pm];
    geom    = permute(geom, pm);
end

% Find all nonzero elements in the geometry array (linear index form).
idxNZ = find(geom);

% Take the differences of the linear indices.
% +1 ==> consecutive node or array-bounded trailing node
% >1 ==> geometry-bounded trailing node
cond = diff(idxNZ);

% Identify geometry-bounded trailing nodes with a 2. the last
% node with a zero as well (forcibly a trailing node).
cond(cond > 1) = 2;

% Identify the last node, which is forcibly a trailing node, with a 2.
cond = [cond; 2];

% Identify array-bounded trailing nodes with a 2.
cond(mod(idxNZ,sz(dir)) == 0) = 2;

% Identify the leading nodes with a 0. The first node is forcibly a leading
% node. Nodes immediately following trailing nodes are also leading nodes.
tmp                        = (cond == 2);
cond([true; tmp(1:end-1)]) = 0;
clear tmp;

% Define the leading, interior and trailing nodes.
idxl = idxNZ(cond == 0);
idxi = idxNZ(cond == 1);
idxt = idxNZ(cond == 2);
clear cond idxNZ;

% Determine the linear indices for the unpermuted geometry array.
if dir > 1
    % Leading node indices.
    subl      = cell(ndims(geom),1);
    [subl{:}] = ind2sub(sz(pm), idxl);
    subl      = [subl(2:pm); subl(1); subl(pm+1:end)];
    idxl      = sub2ind(sz, subl{:});
    clear subl;
    
    % Interior node indices.
    subi      = cell(ndims(geom),1);
    [subi{:}] = ind2sub(sz(pm), idxi);
    subi      = [subi(2:pm); subi(1); subi(pm+1:end)];
    idxi      = sub2ind(sz, subi{:});
    clear subi;
    
    % Trailing node indices.
    subt      = cell(ndims(geom),1);
    [subt{:}] = ind2sub(sz(pm), idxt);
    subt      = [subt(2:pm); subt(1); subt(pm+1:end)];
    idxt      = sub2ind(sz, subt{:});
    clear subt;
end
clear dir pm sz;


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
% 2022/02/25 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
