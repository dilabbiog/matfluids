%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                       DEFORMATION GRADIENT TENSOR                       %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical, Industrial & Aerospace Engineering            %
% Concordia University Montréal, Canada                                   %
%                                                                         %
% Last Update: October 2nd, 2018 by Giuseppe Di Labbio                    %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% SYNTAX                                                                  %
%                                                                         %
% DG = PP_DeformGradTensor(A, ds, scheme);                                %
% DG = PP_DeformGradTensor(A, ds, scheme, mask);                          %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the deformation gradient tensor, a material      %
% displacement gradient tensor, of an advected particle field. The user   %
% may specify the finite difference scheme to be used and whether or not  %
% to use a pre-defined mask.                                              %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'A'          - STRUCT                                                   %
%              - Struct containing fields X and Y holding the positions   %
%                of the advected particles.                               %
% ----------------------------------------------------------------------- %
% 'ds'         - TWO-ELEMENT REAL VECTOR                                  %
%              - Spatial grid spacing in the x and y directions. The user %
%                must specify ds in the form ds = [dx dy] (or             %
%                equivalently ds = [dx; dy]).                             %
% ----------------------------------------------------------------------- %
% 'mask'      - 2D DOUBLE ARRAY                                           %
%             - Mask information for the vector field. Values of 1        %
%               correspond to active nodes while values of 0 correspond   %
%               to blanked nodes.                                         %
% ----------------------------------------------------------------------- %
% 'scheme'    - SPECIFIC STRING                                           %
%             - String denoting which scheme to use. The user can select  %
%               from the explicit second order central scheme ('EXO2'),   %
%               or Lele's second, third or fourth order compact schemes   %
%               ('COMP2', 'COMP3', 'COMP4').                              %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Advect a grid of points placed at each node of a time-dependent double  %
% gyre over the domain (x,y) = [0,2]x[0,1] with a constant grid spacing   %
% of 0.01 over the time interval [0,20] with time-step size 0.1. Use A =  %
% 0.1, epsilon = 0.25, and omega = 2*pi/10. Use the time-step dt = 0.1    %
% for the advection (i.e. use the 'singlestep' option). Compute the       %
% deformation gradient tensor of the advected particle field at the final %
% advection time using an explicit second-order centered finite           %
% difference scheme.                                                      %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> ADV = PP_AdvectGrid(VEC{1}, VEC, t(2), 'singlestep');                %
% >> DG = PP_DeformGradTensor(ADV{21}, [x(2) y(2)], 'EXO2');              %
% >> contourf(VEC{21}.X, VEC{21}.Y, DG.yx, 'EdgeColor', 'None');          %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% FD_COMP2                                                                %
% FD_COMP3                                                                %
% FD_COMP4                                                                %
% FD_EXO2                                                                 %
% FD_SCH                                                                  %
% TDMA                                                                    %
%                                                                         %
% Called in:                                                              %
% PP_CauchyGreenTensor                                                    %
% PP_FTLE                                                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_DeformGradTensor
function [DG] = PP_DeformGradTensor(A, ds, scheme, varargin)

% Determine the number of data points in x (nx) and y (ny).
[ny, nx] = size(A.X);

DG.xx = zeros(ny, nx);
DG.xy = zeros(ny, nx);
DG.yx = zeros(ny, nx);
DG.yy = zeros(ny, nx);

if nargin == 4
    for i = 1:ny
        DG.xx(i,:) = FD_SCH(A.X(i,:), ds(1), scheme, varargin{1}(i,:)).';
        DG.yx(i,:) = FD_SCH(A.Y(i,:), ds(1), scheme, varargin{1}(i,:)).';
    end
    
    for j = 1:nx
        DG.xy(:,j) = FD_SCH(A.X(:,j), ds(2), scheme, varargin{1}(:,j)).';
        DG.yy(:,j) = FD_SCH(A.Y(:,j), ds(2), scheme, varargin{1}(:,j)).';
    end
else
    for i = 1:ny
        DG.xx(i,:) = FD_SCH(A.X(i,:), ds(1), scheme).';
        DG.yx(i,:) = FD_SCH(A.Y(i,:), ds(1), scheme).';
    end
    
    for j = 1:nx
        DG.xy(:,j) = FD_SCH(A.X(:,j), ds(2), scheme).';
        DG.yy(:,j) = FD_SCH(A.Y(:,j), ds(2), scheme).';
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
