%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                        SWIRL STRENGTH (LAMBDA-CI)                       %
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
% Lci = PP_SwirlStrength(VGT);                                            %
% Lci = PP_SwirlStrength(VGT, opt);                                       %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function computes the swirl strength for 2D or 3D data (also known %
% as lambda_ci) as defined by Zhou, Adrian, Balachandar and Kendall [1].  %
% The user can optionally choose to normalize the values of swirl         %
% strength. Recall that the swirl strength is taken as the absolute value %
% of the imaginary part of the complex conjugate pair of eigenvalues      %
% associated with the velocity gradient tensor; therefore the swirl       %
% strength is always a positive quantity and values greater than zero     %
% potentially distinguish regions of local swirling.                      %
%                                                                         %
% References:                                                             %
% [1] Zhou, J., Adrian, R. J., Balachandar, S. & Kendall, T. M. (1999).   %
%     Mechanisms for generating coherent packets of hairpin vortices in   %
%     channel flow. Journal of Fluid Mechanics, 387, 353-396.             %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'Lci'        - 2D DOUBLE ARRAY                                          %
%              - Lambda-ci criterion.                                       %
% ----------------------------------------------------------------------- %
% 'opts'       - SPECIFIC STRING                                          %
%              - Option to normalize the values of lambda_ci using the    %
%                string 'norm'.                                           %
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
% Calculate the swirl strength of a time-dependent double gyre on the     %
% domain (x,y) = [0,2]x[0,1] with a constant grid spacing of 0.01 over    %
% the time interval [0,20] with time-step size 0.1. Use A = 0.1, epsilon  %
% = 0.25, and omega = 2*pi/10. Do not normalize the values of the swirl   %
% strength.                                                               %
%                                                                         %
% >> x = linspace(0, 2, 201).';                                           %
% >> y = linspace(0, 1, 101).';                                           %
% >> t = linspace(0, 20, 21).';                                           %
% >> A = 0.1;                                                             %
% >> epsn = 0.25;                                                         %
% >> omga = 2*pi/10;                                                      %
% >> [VEC, VGT] = GEN_DoubleGyre(x, y, t, A, epsn, omga);                 %
% >> Lci = cell(length(t),1);                                             %
% >> for k = 1:length(t)                                                  %
% Lci{k} = PP_SwirlStrength(VGT{k});                                      %
% end                                                                     %
% >> for k = 1:length(t)                                                  %
% contourf(VEC{k}.X, VEC{k}.Y, Lci{k}, 'EdgeColor', 'None');              %
% colormap gray;                                                          %
% hold on;                                                                %
% quiver(VEC{k}.X(1:4:end,1:4:end), VEC{k}.Y(1:4:end,1:4:end), ...        %
%        VEC{k}.U(1:4:end,1:4:end), VEC{k}.V(1:4:end,1:4:end), 'b');      %
% hold off;                                                               %
% axis([0 2 0 1]);                                                        %
% pause(0.25);                                                            %
% end                                                                     %
% >> clear k;                                                             %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% PP_TensorEig                                                            %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_SwirlStrength
function [Lci] = PP_SwirlStrength(VGT, varargin)

if length(fieldnames(VGT)) == 4
    
    % Solve the tensor eigensystem for the eigenvalues and eigenvectors of
    % the velocity gradient tensor.
    [E1, E2] = PP_TensorEig(VGT, 'none', logical(VGT.UX));
    
    % Compute the swirl strength (lambda_ci) as largest imaginary part of
    % the two eigenvalues. In 2D, the velocity gradient tensor will have
    % either two real eigenvalues or a pair of complex conjugate
    % eigenvalues, with the value of lambda_ci being zero in the real
    % eigenvalue case; here the positive value of the conjugate pairs are
    % taken for the value of lambda_ci. Note also that in 2D, lambda_ci
    % could have been computed using:
    % Lci = 0.5*imag(sqrt((VGT.UX - VGT.VY).^2 + 4*VGT.UY.*VGT.VX));
    Lci = max(imag(E1.eVal), imag(E2.eVal));
    
    % Normalize the swirl strength.
    if nargin == 2 && strcmpi(varargin{1},'norm')
        Lci = Lci/max(max(Lci));
    end
    
elseif length(fieldnames(VGT)) == 9
    
    % Solve the tensor eigensystem for the eigenvalues and eigenvectors of
    % the velocity gradient tensor.
    [E1, E2, E3] = PP_TensorEig(VGT, 'none', logical(VGT.UX));
    
    % Compute the swirl strength (lambda_ci) as largest imaginary part of
    % the three eigenvalues. In 3D, the velocity gradient tensor will have
    % either three real eigenvalues or one real eigenvalue and a pair of
    % complex conjugate eigenvalues, with the value of lambda_ci being zero
    % in the real eigenvalue case; here the positive value of the conjugate
    % pairs are taken for the value of lambda_ci.
    Lci = max(imag(E1.eVal), imag(E2.eVal), imag(E3.eVal));
  
    % Normalize the swirl strength.
    if nargin == 2 && strcmpi(varargin{1},'norm')
        Lci = Lci/max(max(max(Lci)));
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
