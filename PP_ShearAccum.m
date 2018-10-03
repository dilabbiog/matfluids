%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             POST-PROCESSING                             %
%                     SHEAR ACCUMULATION FOR PARTICLES                    %
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
% ShAcc = PP_ShearAccum(P, VEC, VGT, mu, dt);                             %
%                                                                         %
% DESCRIPTION:                                                            %
%                                                                         %
% N/A                                                                     %
%                                                                         %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% N/A                                                                     %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% N/A                                                                     %
%                                                                         %
% DEPENDENCIES                                                            %
%                                                                         %
% Requires:                                                               %
% PP_AdvectPoints                                                         %
% PP_NewtStressTensor                                                     %
%                                                                         %
% Called in:                                                              %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PP_ShearAccum
function [ShAcc] = PP_ShearAccum(P, VEC, VGT, mu, dt)

n = length(VEC);
p = length(P.X);

A = PP_AdvectPoints(P, VEC, dt, 'singlestep');

NST = cell(n,1);
VMS = cell(n,1);
for k = 1:n
    NST{k} = PP_NewtStressTensor(VGT{k},mu);
    VMS{k} = sqrt(NST{k}.xx.^2 - NST{k}.xx.*NST{k}.yy + NST{k}.yy.^2 + ...
                  3*NST{k}.xy.^2)/sqrt(3);
end

ShAcc = zeros(n, p);
for k = 2:n
    ShAcc(k,:) = dt*interp2(VEC{k}.X, VEC{k}.Y, VMS{k}, A.X(k,:), ...
                            A.Y(k,:));
end
ShAcc = cumsum(ShAcc);

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
