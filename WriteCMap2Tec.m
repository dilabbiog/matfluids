%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                      WRITE COLORMAP TO TECPLOT 360                      %
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
% WriteCMap2Tec(CMname);                                                  %
% WriteCMap2Tec(CMname, Name, Value);                                     %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function was sparked by Lukasz Panek's script (publicly available  %
% at [1]) which writes several of MATLAB's colormaps to files for use in  %
% TecPlot 360. My original intention and use for this function was to     %
% translate the perceptually uniform colormaps created by Stéfan van der  %
% Walt and Nathaniel Smith, which were posted on MATLAB's File Exchange   %
% by Ander Biguri at [2], for use in TecPlot 360. I eventually started    %
% making my own equation-governed colormaps in MATLAB which I would bring %
% over to TecPlot using this function.                                    %
%                                                                         %
% Here, any M-by-3 colormap matrix defined in MATLAB can be written to a  %
% '.map' file for use by TecPlot 360. The main benefit of this function   %
% is that the number of control points defined in TecPlot 360 can be      %
% specified by the user for up to 50 points, which is the maximum number  %
% of points allowed by TecPlot. The values at the control points are      %
% interpolated for using MATLAB's 'interp1' function, the method of which %
% can be specified by the user. The user may additionally specify the     %
% name of the saved file and the TecPlot macro version number.            %
%                                                                         %
% References:                                                             %
% [1] http://www.cfd.tu-berlin.de/~panek/cfd/colormaps/teccolormap.m      %
% [2] https://www.mathworks.com/matlabcentral/fileexchange/51986-         %
%     perceptually-uniform-colormaps                                      %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'CMname'     - Name of the colormap to export from MATLAB to TecPlot.   %
%              - STRING                                                   %
% ----------------------------------------------------------------------- %
% Name-Value Pairs:                                                       %
% ----------------------------------------------------------------------- %
% 'FileName'   - Name of the exported colormap file.                      %
%              - STRING                                                   %
%              - Default: 'MATLABcmap'                                    %
% ----------------------------------------------------------------------- %
% 'MacroVer'   - Macro version of TecPlot 360.                            %
%              - FOUR DIGIT STRING                                        %
%              - Default: '1410'                                          %
% ----------------------------------------------------------------------- %
% 'Method'     - Method of interpolation to use (for interp1).            %
%              - 'linear', 'nearest', 'next', 'previous', 'pchip',        %
%                'cubic' , 'spline'                                       %
%              - Default: 'cubic'                                         %
% ----------------------------------------------------------------------- %
% 'NumCtrlPts' - Number of control points to use for TecPlot colormap.    %
%              - INTEGER (0 < N <= 50)                                    %
%              - Default: 20                                              %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE 1                                                               %
%                                                                         %
% Export the perceptually uniform colormap 'viridis' for TecPlot 360      %
% using 30 control points and linear interpolation.                       %
% >> WriteCMap2Tec('viridis', 'NumCtrlPts', 30, 'Method', 'linear');      %
%                                                                         %
% EXAMPLE 2                                                               %
%                                                                         %
% Export the 'parula' colormap for TecPlot 360 using 15 control points    %
% and name the file 'MATparula'.                                          %
% >> WriteCMap2Tec('parula', 'NumCtrlPts', 15, 'FileName', 'MATparula');  %
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

%% WriteCMap2Tec
function [] = WriteCMap2Tec(CMname, varargin)

%% Parse the inputs.
chk = {'filename', 'macrover', 'method', 'numctrlpts'};
hasChoice  = struct('fname', 0, 'mver', 0, 'method', 0, 'ctrlpts', 0);
numCalls   = [0 0 0 0];
methodOpts = {'linear', 'nearest', 'next', 'previous', 'pchip',         ...
              'cubic', 'spline'};
if (nargin >= 1) && (nargin <= 9) && mod(nargin,2)
    if nargin > 1
        for k = 1:2:nargin-1
            switch lower(varargin{k})
                case chk{1}
                    cond =  ischar(varargin{k+1})                       ...
                         & ~max(strcmpi(chk, varargin{k+1}));
                    if cond
                        hasChoice.fname   = k + 1;
                        numCalls(1)       = numCalls(1) + 1;
                    else
                        error('InputError:Fname',                       ...
                              'The file name is not valid');
                    end
                case chk{2}
                    cond =  ischar(varargin{k+1})                       ...
                         & (length(varargin{k+1}) == 4)                 ...
                         &  min(isstrprop(varargin{k+1}, 'digit'));
                    if cond
                        hasChoice.mver    = k + 1;
                        numCalls(2)       = numCalls(2) + 1;
                    else
                        error('InputError:Mver',                        ...
                              'The macro version is not valid');
                    end
                case chk{3}
                    cond =  ischar(varargin{k+1})                       ...
                         &  max(strcmpi(methodOpts, varargin{k+1}));
                    if cond
                        hasChoice.method  = k + 1;
                        numCalls(3)       = numCalls(3) + 1;
                    else
                        error('InputError:Method',                      ...
                              'The selected method is not valid');
                    end
                case chk{4}
                    cond = (max(size(varargin{k+1})) == 1)              ...
                         &  isnumeric(varargin{k+1})                    ...
                         &  isreal(varargin{k+1})                       ...
                         & (varargin{k+1} > 0)                          ...
                         & (varargin{k+1} <= 50)                        ...
                         & (floor(varargin{k+1}) == varargin{k+1});
                    if cond
                        hasChoice.ctrlpts = k + 1;
                        numCalls(4)       = numCalls(4) + 1;
                    else
                        error('InputError:CtrlPts',                     ...
                              'The number of control points is not valid');
                    end
            end
            if ~max(strcmpi(varargin{k}, chk))
                error('InputError:InvalidInputs',                       ...
                      'One or more inputs are not recognized.');
            end
            if max(numCalls) > 1
                error('InputError:Repeated',                            ...
                      'One or more inputs have been repeated.');
            end
        end
    end
else
    error('InputError:NumOpts', 'Invalid number of input arguments.');
end
clear chk cond k methodOpts numCalls;

%% Rename the inputs for ease of access.
if hasChoice.fname  , fname  = varargin{hasChoice.fname  };
else                  fname  = 'MATLABcmap';                end
if hasChoice.mver   , mver   = varargin{hasChoice.mver   };
else                  mver   = '1410';                      end
if hasChoice.method , method = varargin{hasChoice.method };
else                  method = 'cubic';                     end
if hasChoice.ctrlpts, Ncp    = varargin{hasChoice.ctrlpts};
else                  Ncp    = 20;                          end

%% Load the colormap.
map = colormap(CMname);
close all;

%% Determine the number points in the colormap.
pts = size(map, 1);
if Ncp > pts
    error('InputError:TooManyPts',                                      ...
         ['More control points have been selected than are available '  ...
          'in the colormap.']);
end

%% Determine the colormap values at the control points.
OldPts = (0:pts-1).';
NewPts = linspace(0, pts-1, Ncp).';
NewMap = zeros(Ncp, 3);

NewMap(:,1) = interp1(OldPts, map(:,1), NewPts, lower(method));
NewMap(:,2) = interp1(OldPts, map(:,2), NewPts, lower(method));
NewMap(:,3) = interp1(OldPts, map(:,3), NewPts, lower(method));
clear map method OldPts;

%% Open a file for writing the colormap.
fid = fopen(sprintf('%s.map', fname), 'w');
if fid == -1
    error('OutputFile:CantOpen', 'Can''t open colormap file for output.')
end
clear fname;

%% Write the file header.
fprintf(fid, '#!MC %s\r\n', mver);
fprintf(fid, '\r\n');
fprintf(fid, '$!COLORMAP\r\n');
fprintf(fid, '  CONTOURCOLORMAP = USERDEF\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '$!COLORMAPCONTROL RESETTOFACTORY\r\n');
fprintf(fid, '$!COLORMAP\r\n');
fprintf(fid, '  USERDEFINED\r\n');
fprintf(fid, '  {\r\n');
fprintf(fid, '   NUMCONTROLPOINTS = %d\r\n', Ncp);
clear mver;

%% Write the control points to the file.
for k = 1:Ncp
    fprintf(fid, '     CONTROLPOINT %d\r\n', k);
    fprintf(fid, '     {\r\n');
    fprintf(fid, '      COLORMAPFRACTION = %1.6f\r\n', NewPts(k)/(pts-1));
    fprintf(fid, '      LEADRGB\r\n');
    fprintf(fid, '      {\r\n');
    fprintf(fid, '       R = %d\r\n', round(NewMap(k,1)*255));
    fprintf(fid, '       G = %d\r\n', round(NewMap(k,2)*255));
    fprintf(fid, '       B = %d\r\n', round(NewMap(k,3)*255));
    fprintf(fid, '      }\r\n');
    fprintf(fid, '      TRAILRGB\r\n');
    fprintf(fid, '      {\r\n');
    fprintf(fid, '       R = %d\r\n', round(NewMap(k,1)*255));
    fprintf(fid, '       G = %d\r\n', round(NewMap(k,2)*255));
    fprintf(fid, '       B = %d\r\n', round(NewMap(k,3)*255));
    fprintf(fid, '      }\r\n');
    fprintf(fid, '     }\r\n');
end
fprintf(fid, '  }\r\n');
clear k Ncp NewMap NewPts pts;

%% Close the file.
fclose(fid);
clear fid;

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*N/A>
% Line(s) N/A
% Message(s)
% * N/A.
% Reason(s)
% * N/A.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
