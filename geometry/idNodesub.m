%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            GEOMETRY TOOLBOX                             %
%                                                                         %
% idNodesub                                                               %
% Node Property Identification                                            %
% Examine geometry for subscripts of leading, trailing and interior nodes %
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
% [subl, subi, subt] = idNodesub(geom);                                   %
% [subl, subi, subt] = idNodesub(geom, drc);                              %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Determine the subscripts of the leading, interior and trailing nodes of %
% an input array. By default, the subscripts are listed in the order they %
% appear along the first dimension (since this is how linear indices are  %
% counted by MATLAB). Alternatively, the subscripts may be listed in the  %
% order they appear along a different direction.                          %
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
% 'drc'          POSITIVE INTEGER SCALAR                                  %
%                Default: 1                                               %
%              ~ Input scalar. Direction along which to list the          %
%                subscripts.                                              %
% ======================================================================= %
% Name-Value Pair Arguments:                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ======================================================================= %
% Output Arguments:                                                       %
% ----------------------------------------------------------------------- %
% 'subl'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of leading nodes in the %
%                order they appear along a specified direction.           %
% ----------------------------------------------------------------------- %
% 'subi'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of interior nodes in    %
%                the order they appear along a specified direction.       %
% ----------------------------------------------------------------------- %
% 'subt'         CELL ARRAY                                               %
%              ~ Output cell array of size ndims(geom) by 1. Each cell k  %
%                contains the list of subscript k of trailing nodes in    %
%                the order they appear along a specified direction.       %
% ======================================================================= %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in one dimension.                                               %
%                                                                         %
% >> geom               = [ones(20,1); zeros(10,1); ones(6,1)];           %
% >> [subl, subi, subt] = idNodesub(geom);                                %
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
% >> [subl, subi, subt] = idNodesub(geom, 1);                             %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                         %
% EXAMPLE 3                                                               %
%                                                                         %
% Determine the leading, interior and trailing nodes for a geometry       %
% defined in two dimensions. List the subscripts in the order they appear %
% along the second dimension.                                             %
%                                                                         %
% >> geom               = zeros(51, 101);                                 %
% >> geom(5:20,5:20)    = 1;                                              %
% >> geom(40:45,55:80)  = 1;                                              %
% >> [subl, subi, subt] = idNodesub(geom, 2);                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [subl, subi, subt] = idNodesub(geom, varargin)


%% PARSE INPUTS

% Input defaults.
default.drc = 1;

% Input checks.
check.drc  = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'finite', 'positive', 'scalar'});
check.geom = @(x) validateattributes(x,                                 ...
                  {'logical', 'numeric'},                               ...
                  {'nonempty', 'binary'});

% Parse the inputs.
hParser = inputParser;
addRequired ( hParser, 'geom' , check.geom              );
addOptional ( hParser, 'drc'  , default.drc , check.drc );
parse(hParser, geom, varargin{:});
clear check default;

% Additional verifications.
nargoutchk(0,3);


%% INDEX SEARCH

% Determine the size and number of dimensions of the input array.
sz = size(geom);
nd = ndims(geom);

% Define the search direction and its length.
drc = hParser.Results.drc;
len = sz(drc);

% Permute the input array.
pm = 1:nd;
if drc > 1
    pm(drc) = [];
    pm      = [drc pm];
    geom    = permute(geom, pm);
end
szp = sz(pm);

% Find all nonzero elements in the input array (linear index form).
idxNZ = find(geom);

% Take the differences of the linear indices.
% +1 ==> consecutive node or array-bounded trailing node
% >1 ==> geometry- or array-bounded trailing node
id = diff(idxNZ);

% Identify known trailing nodes with a 2.
id(id > 1) = 2;

% Identify the last node, which is forcibly a trailing node, with a 2.
id = [id; 2];

% Identify array-bounded trailing nodes with a 2.
id(mod(idxNZ,len) == 0) = 2;

% Identify the leading nodes with a 0. The first node is forcibly a leading
% node. Nodes immediately following trailing nodes are also leading nodes.
tmp                      = (id == 2);
id([true; tmp(1:end-1)]) = 0;
clear tmp;

% Define the leading, interior and trailing nodes.
idxl = idxNZ(id == 0);
idxi = idxNZ(id == 1);
idxt = idxNZ(id == 2);
clear id idxNZ;

% Determine the leading node subscripts.
subl      = cell(nd,1);
[subl{:}] = ind2sub(szp, idxl);
clear idxl;

% Determine the interior node subscripts.
subi      = cell(nd,1);
[subi{:}] = ind2sub(szp, idxi);
clear idxi;

% Determine the trailing node subscripts.
subt      = cell(nd,1);
[subt{:}] = ind2sub(szp, idxt);
clear idxt;

% Determine the subscripts for the unpermuted input array.
if drc > 1
    subl = [subl(2:pm); subl(1); subl(pm+1:end)];
    subi = [subi(2:pm); subi(1); subi(pm+1:end)];
    subt = [subt(2:pm); subt(1); subt(pm+1:end)];
end
clear drc len nd pm sz szp;


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
% 2022/02/28 -- (GDL) Beta version of the code finalized.                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% FUTURE UPDATES                                                          %
%                                                                         %
% None foreseen at the moment.                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
