%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                               Q-CRITERION                               %
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
% SYNTAX                                                                  %
%                                                                         %
% Q = PP_InvariantQ(VGT);                                                 %
% Q = PP_InvariantQ(VGT, opts);                                           %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the Q-criterion, the second invariant of the     %
% velocity gradient tensor, for 2D or 3D data as defined by Hunt, Wray    %
% Moin [1]. The user can optionally choose to keep only positive values   %
% of Q and to normalize the values of Q with respect to its largest       %
% value. Recall that positive values of Q signify the presence of a       %
% vortex.                                                                 %
%                                                                         %
% References:                                                             %
% [1] Hunt, J. C. R., Wray, A. A. & Moin, P. (1988). Eddies, streams, and %
%     convergence zones in turbulent flows. Proceedings of the Summer     %
%     Program of the Center for Turbulent Research, 193-208.              %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'Q'          - 2D DOUBLE ARRAY                                          %
%              - Q criterion (second invariant of the velocity gradient   %
%                tensor).                                                 %
% ----------------------------------------------------------------------- %
% 'opts'       - SPECIFIC STRING                                          %
%              - Option to normalize the values of Q by its maximum value %
%                (not its absolute maximum) using 'norm', to retain only  %
%                positive values of Q using 'pos', or both.               %
% ----------------------------------------------------------------------- %
% 'VGT'        - STRUCT                                                   %
%              - Velocity gradient tensor containing the four components  %
%                'UX', 'UY', 'VX' and 'VY' in 2D or the nine components   %
%                'UX', 'UY', 'UZ', 'VX', 'VY', 'VZ', 'WX', 'WY' and 'WZ'  %
%                in 3D.                                                   %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Calculate the Q-criterion of a time-dependent double gyre on the domain %
% (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over the time  %
% interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon = 0.25,   %
% and omega = 2*pi/10. Keep only positive values of normalized Q.         %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> Q = cell(length(t),1);                                               %
% >> for k = 1:length(t)                                                  %
% Q{k} = PP_InvariantQ(VGT{k}, 'pos', 'norm');                            %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, Q{k}, 'EdgeColor', 'None');                %
% colormap gray;                                                          %
% hold on;                                                                %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'w');      %
% hold off;                                                               %
% axis([0 2 0 1]);                                                        %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% N/A                                                                     %
%                                                                         %
% Called in:                                                              %
% PP_DeltaCrit                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_InvariantQ
function [Q] = PP_InvariantQ(VGT, varargin)

if length(fieldnames(VGT)) == 4
    
    % Compute the Q criterion.
    Q = -0.5*(VGT.UX.^2 + 2*VGT.UY.*VGT.VX + VGT.VY.^2);
    
    % Apply the normalization and/or keep the positive values only.
    switch nargin
        case 2
            if strcmpi(varargin{1},'pos')
                Q(Q < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                Q = Q/max(max(Q));
            end
        case 3
            if strcmpi(varargin{1},'pos') || strcmpi(varargin{2},'pos')
                Q(Q < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                Q = Q/max(max(Q));
            end
    end
    
elseif length(fieldnames(VGT)) == 9
    
    % Compute the Q criterion.
    Q = -0.5*(VGT.UX.^2 + 2*VGT.UY.*VGT.VX + 2*VGT.UZ.*VGT.WX ...
      +       VGT.VY.^2 + 2*VGT.VZ.*VGT.WY ...
      +       VGT.WZ.^2);
    
    % Apply the normalization and/or keep the positive values only.
    switch nargin
        case 2
            if strcmpi(varargin{1},'pos')
                Q(Q < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                Q = Q/max(max(max(Q)));
            end
        case 3
            if strcmpi(varargin{1},'pos') || strcmpi(varargin{2},'pos')
                Q(Q < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                Q = Q/max(max(max(Q)));
            end
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
