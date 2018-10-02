%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                             DELTA-CRITERION                             %
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
% D = PP_DeltaCrit(VGT);                                                  %
% D = PP_DeltaCrit(VGT, opts);                                            %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the delta criterion for 2D or 3D data as defined %
% by Chong, Perry and Cantwell (see Chong, M. S., Perry, A. E. &          %
% Cantwell, B. J. (1990). A general classification of three-dimensional   %
% flow fields. Physics of Fluids A: Fluid Dynamics, 2, 765-777). The user %
% can optionally choose to keep only positive values of delta and to      %
% normalize the values of delta with respect to its largest value. Recall %
% that positive values of delta signify the presence of a vortex.         %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'D'          - 2D DOUBLE ARRAY                                          %
%              - Delta criterion.                                         %
% ----------------------------------------------------------------------- %
% 'opts'       - SPECIFIC STRING                                          %
%              - Option to normalize the values of delta by its maximum   %
%                value (not its absolute maximum) using 'norm', to retain %
%                only positive values of delta using 'pos', or both.      %
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
% Calculate the delta criterion of a time-dependent double gyre on the    %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10. Keep only positive values of normalized    %
% delta.                                                                  %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> D = cell(length(t),1);                                               %
% >> for k = 1:length(t)                                                  %
% D{k} = PP_DeltaCrit(VGT{k}, 'pos', 'norm');                             %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, D{k}, 'EdgeColor', 'None');                %
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
% PP_InvariantQ                                                           %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_DeltaCrit
function [D] = PP_DeltaCrit(VGT, varargin)

% Compute the second invariant of the velocity gradient tensor.
Q = PP_InvariantQ(VGT);

if length(fieldnames(VGT)) == 4
    
    % Determine the size of the flow domain.
    [ny,nx] = size(VGT.UX);
    
    % Compute the determinant of the velocity gradient tensor.
    detVGT = zeros(ny,nx);
    for i = 1:ny
        for j = 1:nx
            detVGT(i,j) = det([VGT.UX(i,j) VGT.UY(i,j);
                               VGT.VX(i,j) VGT.VY(i,j)]);
        end
    end
    
    % Compute the delta criterion.
    D = (Q/3).^3 + (detVGT/2).^2;
    
    % Apply the normalization and/or keep the positive values only.
    switch nargin
        case 2
            if strcmpi(varargin{1},'pos')
                D(D < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                D = D/max(max(D));
            end
        case 3
            if strcmpi(varargin{1},'pos') || strcmpi(varargin{2},'pos')
                D(D < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                D = D/max(max(D));
            end
    end
    
elseif length(fieldnames(VGT)) == 9
    
    % Determine the size of the flow domain.
    [ny,nx,nz] = size(VGT.UX);
    
    % Compute the determinant of the velocity gradient tensor.
    detVGT = zeros(ny,nx,nz);
    for i = 1:ny
        for j = 1:nx
            for k = 1:nz
                detVGT(i,j,k) = det([VGT.UX(i,j,k)  VGT.UY(i,j,k)
                                     VGT.UZ(i,j,k); VGT.VX(i,j,k)
                                     VGT.VY(i,j,k)  VGT.VZ(i,j,k);
                                     VGT.WX(i,j,k)  VGT.WY(i,j,k)
                                     VGT.WZ(i,j,k)]);
            end
        end
    end
    
    % Compute the delta criterion.
    D = (Q/3).^3 + (detVGT/2).^2;
    
    % Apply the normalization and/or keep the positive values only.
    switch nargin
        case 2
            if strcmpi(varargin{1},'pos')
                D(D < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm')
                D = D/max(max(max(D)));
            end
        case 3
            if strcmpi(varargin{1},'pos') || strcmpi(varargin{2},'pos')
                D(D < 0) = 0;
            end
            
            if strcmpi(varargin{1},'norm') || strcmpi(varargin{2},'norm')
                D = D/max(max(max(D)));
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
