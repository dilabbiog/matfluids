%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                        WRITE DATA TO TECPLOT 360                        %
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
% Write2TecPlot(fname, T, varnames, VEC, ...);                            %
% Write2TecPlot(fname, T, varnames, X, Y, ...);                           %
% Write2TecPlot(fname, T, varnames, X, Y, U, V, ...);                     %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function writes data to files that can be read by TecPlot 360.     %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'fname'      - STRING                                                   %
%              - File name of written data.                               %
% ----------------------------------------------------------------------- %
% 'T'          - 1D DOUBLE ARRAY                                          %
%              - Time vector of the velocity field.                       %
% ----------------------------------------------------------------------- %
% 'U'          - 1D CELL, ELEMENTS: 2D DOUBLE ARRAYS                      %
%              - Component of velocity vector field in the x direction in %
%                the form of a cell array.                                %
% ----------------------------------------------------------------------- %
% 'V'          - 1D CELL, ELEMENTS: 2D DOUBLE ARRAYS                      %
%              - Component of velocity vector field in the y direction in %
%                the form of a cell array.                                %
% ----------------------------------------------------------------------- %
% 'varnames'   - 1D CELL, ELEMENTS: STRINGS                               %
%              - List of names of additional variables in the form of a   %
%                cell array of strings. If no additional variables will   %
%                written, use varnames = {}.                              %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS, ELEMENTS: 2D DOUBLE ARRAYS   %
%              - Velocity field data containing the mask 'C', Cartesian   %
%                coordinates 'X', 'Y' and 'Z', and velocity field         %
%                components 'U', 'V' and 'W'.                             %
% ----------------------------------------------------------------------- %
% 'X'          - 1D CELL, ELEMENTS: 2D DOUBLE ARRAYS                      %
%              - Gridded Cartesian coordinate x in the form of a cell     %
%                array.                                                   %
% ----------------------------------------------------------------------- %
% 'Y'          - 1D CELL, ELEMENTS: 2D DOUBLE ARRAYS                      %
%              - Gridded Cartesian coordinate y in the form of a cell     %
%                array.                                                   %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Write the VEC variable, read from DaVis data, to TecPlot 360.           %
%                                                                         %
% >> Write2TecPlot('Flow', T, {}, VEC);                                   %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Write the VEC variable, read from DaVis data, as well as a computed     %
% vorticity field to TecPlot 360.                                         %
%                                                                         %
% >> n = length(VEC);                                                     %
% >> for k = 1:n                                                          %
% VORT{k} = PP_Vorticity(VEC{k});                                         %
% end                                                                     %
% >> Write2TecPlot('Flow', T, {'vort'}, VEC, VORT);                       %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% mat2vec                                                                 %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write2TecPlot
function [] = Write2TecPlot(fname, T, varnames, varargin)

for k = 4:nargin
    % Check if the variable VEC is one of the input arguments.
    if strcmpi('VEC', inputname(k)) 
        hasVEC = k-3;
        break;
    elseif isstruct(varargin{k-3}{1}) && (isfield(varargin{k-3}{1},'X') ...
                                      &&  isfield(varargin{k-3}{1},'Y') ...
                                      &&  isfield(varargin{k-3}{1},'U') ...
                                      &&  isfield(varargin{k-3}{1},'V'))
        hasVEC = k-3;
        break;
    else
        hasVEC = 0;
    end
end

hasX = 0;
hasY = 0;
hasU = 0;
hasV = 0;
for k = 4:nargin
    % Check if the variable X is one of the input arguments.
    if ~hasX && strcmpi('X', inputname(k))
        hasX = k-3;
    end
    
    % Check if the variable Y is one of the input arguments.
    if ~hasY && strcmpi('Y', inputname(k))
        hasY = k-3;
    end
    
    % Check if the variable U is one of the input arguments.
    if ~hasU && strcmpi('U', inputname(k))
        hasU = k-3;
    end
    
    % Check if the variable V is one of the input arguments.
    if ~hasV && strcmpi('V', inputname(k))
        hasV = k-3;
    end
    
end

if hasVEC && (~hasX && ~hasY && ~hasU && ~hasV)
    % Determine the number of time steps in the time series.
    n = length(varargin{hasVEC});
    
    % Determine the number of additional arguments.
    nvars = nargin - 4;
    
    for k = 1:n
        
        % Determine the domain size of the velocity field.
        [ny,nx] = size(varargin{hasVEC}{k}.('X'));
        
        % Convert the x, y, u and v maps into lists.
        x = mat2vec(varargin{hasVEC}{k}.('X').');
        y = mat2vec(varargin{hasVEC}{k}.('Y').');
        u = mat2vec(varargin{hasVEC}{k}.('U').');
        v = mat2vec(varargin{hasVEC}{k}.('V').');
        
        % Open a new file.
        fid = fopen(sprintf([fname '_%05d.dat'], k), 'w');
        
        % Type the file header for TecPlot recognition.
        fprintf(fid, 'TITLE = "MATLAB Exported Data"\r\n');
        fprintf(fid, 'FILETYPE = FULL\r\n');
        fprintf(fid, 'VARIABLES = "X", "Y", "U", "V"');
        if nvars
            for p = 1:nvars
                fprintf(fid, ', "%s"', varnames{p});
            end
        end
        fprintf(fid, '\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.UVar = "3"\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.VVar = "4"\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.VectorVarsAreVelocity = "TRUE"\r\n');
        fprintf(fid, 'ZONE\r\n');
        fprintf(fid, 'T = "TimeStep_%05d"\r\n', k);
        fprintf(fid, 'STRANDID = 1, SOLUTIONTIME = %5.3f\r\n', T(k));
        fprintf(fid, 'I = %d, J = %d, K = 1\r\n', nx, ny);
        fprintf(fid, 'ZONETYPE = Ordered, DATAPACKING = POINT\r\n');
        fprintf(fid, 'DT = (SINGLE SINGLE SINGLE SINGLE');
        if nvars
            for p = 1:nvars
                fprintf(fid, ' SINGLE');
            end
        end
        fprintf(fid, ')\r\n');
        data = [x, y, u, v];
        if nvars
            for p = 1:nvars+1
                if p == hasVEC, continue; end
                data = [data, mat2vec(varargin{p}{k}.')];
            end
        end
        fprintf(fid, ['%12.8f %12.8f %12.8f %12.8f' ...
                      repmat(' %12.8f', 1, nvars) '\r\n'], data.');
        fprintf(fid, '\r\n');
        
        % Close the file.
        fclose(fid);
        
    end
    
elseif (hasX && hasY) && (~hasVEC && ~hasU && ~hasV)
    % Determine the number of time steps in the time series.
    n = length(varargin{hasX});
    
    % Determine the number of additional arguments.
    nvars = nargin - 5;
    
    for k = 1:n
        
        % Determine the domain size of the velocity field.
        [ny,nx] = size(varargin{hasX}{k});
        
        % Convert the x and y maps into lists.
        x = mat2vec(varargin{hasX}{k}.');
        y = mat2vec(varargin{hasY}{k}.');
        
        % Open a new file.
        fid = fopen(sprintf([fname '_%05d.dat'], k), 'w');
        
        % Type the file header for TecPlot recognition.
        fprintf(fid, 'TITLE = "MATLAB Exported Data"\r\n');
        fprintf(fid, 'FILETYPE = FULL\r\n');
        fprintf(fid, 'VARIABLES = "X", "Y"');
        if nvars
            for p = 1:nvars
                fprintf(fid, ', "%s"', varnames{p});
            end
        end
        fprintf(fid, '\r\n');
        fprintf(fid, 'ZONE\r\n');
        fprintf(fid, 'T = "TimeStep_%05d"\r\n', k);
        fprintf(fid, 'STRANDID = 1, SOLUTIONTIME = %5.3f\r\n', T(k));
        fprintf(fid, 'I = %d, J = %d, K = 1\r\n', nx, ny);
        fprintf(fid, 'ZONETYPE = Ordered, DATAPACKING = POINT\r\n');
        fprintf(fid, 'DT = (SINGLE SINGLE');
        if nvars
            for p = 1:nvars
                fprintf(fid, ' SINGLE');
            end
        end
        fprintf(fid, ')\r\n');
        data = [x, y];
        if nvars
            for p = 1:nvars+2
                if (p == hasX || p == hasY), continue; end
                data = [data, mat2vec(varargin{p}{k}.')];
            end
        end
        fprintf(fid, ['%12.8f %12.8f' ...
                      repmat(' %12.8f', 1, nvars) '\r\n'], data.');
        fprintf(fid, '\r\n');
        
        % Close the file.
        fclose(fid);
        
    end
    
elseif (hasX && hasY && hasU && hasV) && ~hasVEC
    % Determine the number of time steps in the time series.
    n = length(varargin{hasX});
    
    % Determine the number of additional arguments.
    nvars = nargin - 7;
    
    for k = 1:n
        
        % Determine the domain size of the velocity field.
        [ny,nx] = size(varargin{hasX}{k});
        
        % Convert the x, y, u and v maps into lists.
        x = mat2vec(varargin{hasX}{k}.');
        y = mat2vec(varargin{hasY}{k}.');
        u = mat2vec(varargin{hasU}{k}.');
        v = mat2vec(varargin{hasV}{k}.');
        
        % Open a new file.
        fid = fopen(sprintf([fname '_%05d.dat'], k), 'w');
        
        % Type the file header for TecPlot recognition.
        fprintf(fid, 'TITLE = "MATLAB Exported Data"\r\n');
        fprintf(fid, 'FILETYPE = FULL\r\n');
        fprintf(fid, 'VARIABLES = "X", "Y", "U", "V"');
        if nvars
            for p = 1:nvars
                fprintf(fid, ', "%s"', varnames{p});
            end
        end
        fprintf(fid, '\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.UVar = "3"\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.VVar = "4"\r\n');
        fprintf(fid, 'DATASETAUXDATA Common.VectorVarsAreVelocity = "TRUE"\r\n');
        fprintf(fid, 'ZONE\r\n');
        fprintf(fid, 'T = "TimeStep_%05d"\r\n', k);
        fprintf(fid, 'STRANDID = 1, SOLUTIONTIME = %5.3f\r\n', T(k));
        fprintf(fid, 'I = %d, J = %d, K = 1\r\n', nx, ny);
        fprintf(fid, 'ZONETYPE = Ordered, DATAPACKING = POINT\r\n');
        fprintf(fid, 'DT = (SINGLE SINGLE');
        if nvars
            for p = 1:nvars
                fprintf(fid, ' SINGLE');
            end
        end
        fprintf(fid, ')\r\n');
        data = [x, y, u, v];
        if nvars
            for p = 1:nvars+4
                if p == hasX || p == hasY || p == hasU || p == hasV
                    continue;
                end
                data = [data, mat2vec(varargin{p}{k}.')];
            end
        end
        fprintf(fid, ['%12.8f %12.8f %12.8f %12.8f' ...
                      repmat(' %12.8f', 1, nvars) '\r\n'], data.');
        fprintf(fid, '\r\n');
        
        % Close the file.
        fclose(fid);
        
    end
    
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
