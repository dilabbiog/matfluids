%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%          FLOW RECONSTRUCTION BY PROPER ORTHOGONAL DECOMPOSITION         %
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
% VRC = PP_PODreconstruct(VEC, N);                                        %
% VRC = PP_PODreconstruct(VEC, set);                                      %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the proper orthogonal decomposition (POD) of     %
% velocity field data using the algorithm described by L. Sirovich [1].   %
% The algorithm is also well-described by K. Meyer [2]. The user has the  %
% option of performing the POD in mean deviation form where the mean flow %
% is subtracted off at each time step. This function also computes the    %
% global entropy of the flow as defined by N. Aubry [3]. The user can     %
% specify the number of modes to use to reconstruct the data set or       %
% specify exactly which modes to use in the reconstruction.               %
%                                                                         %
% References:                                                             %
% [1] Sirovich, L. (1987). Turbulence and the dynamics of coherent        %
%     structures, Part I: Coherent structures. Quarterly of Applied       %
%     Mathematics, 45(3), 561-571.                                        %
% [2] Meyer, K., Pedersen, J. M. & Özkan, O. (2007). A turbulent jet in   %
%     crossflow analysed with proper orthogonal decomposition. Journal of %
%     Fluid Mechanics, 583, 199-227.                                      %
% [3] Aubry, N. (1991). On the hidden beauty of proper orthogonal         %
%     decomposition. Theoretical and Computational Fluid Dynamics, 2,     %
%     339-352.                                                            %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'N'          - INTEGER SCALAR                                           %
%              - Number of proper orthogonal modes to use in              %
%                reconstructing the data set (uses the first N modes).    %
% ----------------------------------------------------------------------- %
% 'set'        - INTEGER VECTOR                                           %
%              - Integer array specifying exactly which modes to use in   %
%                reconstructing the data set.                             %
% ----------------------------------------------------------------------- %
% 'VEC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Velocity field data containing the mask 'C', Cartesian   %
%                coordinates 'X', 'Y' and 'Z', and velocity field         %
%                components 'U', 'V' and 'W'.                             %
% ----------------------------------------------------------------------- %
% 'VRC'        - 1D CELL, ELEMENTS: STRUCTS                               %
%              - Reconstructed velocity field data containing the mask    %
%                'C', Cartesian coordinates 'X', 'Y' and 'Z', and         %
%                velocity field components 'U', 'V' and 'W'.              %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Reconstruct the velocity field using the first three modes of the       %
% proper orthogonal decomposition of a time-dependent double gyre on the  %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10.                                            %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> VRC = PP_PODreconstruct(VEC, 3);                                     %
% >> for k = 1:length(t)                                                  %
% quiver(VRC{k}.X(1:4:end,1:4:end), VRC{k}.Y(1:4:end,1:4:end), ...        %
%        VRC{k}.U(1:4:end,1:4:end), VRC{k}.V(1:4:end,1:4:end), 'k');      %
% axis([0 2 0 1]);                                                        %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% field2cell                                                              %
% field2mat                                                               %
% PP_POD                                                                  %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_PODreconstruct
function [VRC] = PP_PODreconstruct(VEC, N)

% Compute the proper orthogonal decomposition of the data set.
[POD, ~] = PP_POD(VEC);

% Extract the modal velocity components (U and V) into their own cell
% arrays. Extract also the mode amplitudes and arrange them in a matrix.
U = field2cell(POD, 'U');
V = field2cell(POD, 'V');
a = field2mat(POD, 'a');


% Initialize the reconstructed velocity field.
VRC = cell(length(VEC), 1);
for k = 1:length(VEC)
    VRC{k} = struct('C', VEC{k}.C, 'X', VEC{k}.X, 'Y', VEC{k}.Y, ...
                    'U', 0, 'V', 0);
end

if length(N) == 1
    
    % Reconstruct the velocity field using the first N modes.
    for k = 1:length(VEC)
        for m = 1:N
            
            VRC{k}.U = VRC{k}.U + a(k,m)*U{m};
            VRC{k}.V = VRC{k}.V + a(k,m)*V{m};
            
        end
    end
    
else
    
    % Reconstruct the velocity fields using the specified modes.
    for k = 1:length(VEC)
        for m = 1:length(N)
            
            VRC{k}.U = VRC{k}.U + a(k,N(m))*U{N(m)};
            VRC{k}.V = VRC{k}.V + a(k,N(m))*V{N(m)};
            
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
