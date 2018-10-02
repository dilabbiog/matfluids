%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                         EXTRACT VELOCITY PROFILE                        %
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
% VP = PP_ExtractVelProf(VEC, r);                                         %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% Coming soon ...                                                         %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'N/A'        - N/A                                                      %
%              - N/A                                                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Coming soon ...                                                         %
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

%% PP_ExtractVelProf
function [VP] = PP_ExtractVelProf(VEC, p, r, ds, phi)

VP = struct('X', 0, 'Y', 0, 'U', 0, 'V', 0);

tmp = (-r:ds:r).';
VP.X = tmp*cos(phi) + p(1);
VP.Y = tmp*sin(phi) + p(2);

v = [VP.X(end)-VP.X(1) VP.Y(end)-VP.Y(1)];
v = v/norm(v);
n = [v(2) -v(1)];

Us = interp2(VEC.X, VEC.Y, VEC.U, VP.X, VP.Y, 'cubic');
Vs = interp2(VEC.X, VEC.Y, VEC.V, VP.X, VP.Y, 'cubic');

VP.U = zeros(size(VP.X,1),1);
VP.V = zeros(size(VP.X,1),1);
for k = 1:size(VP.X,1)
    tmp = dot([Us(k) Vs(k)],n)*n;
    VP.U(k) = tmp(1);
    VP.V(k) = tmp(2);
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
