%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GENERAL TOOLBOX                             %
%                    CHECK IF INPUT IS A FORM OF 'YES'                    %
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
% val = isyes(word);                                                      %
%                                                                         %
% DESCRIPTION                                                             %
%                                                                         %
% This function checks whether some input is a form of 'yes'.             %
%                                                                         %
% References:                                                             %
% N/A                                                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
% Variables:                                                              %
% ----------------------------------------------------------------------- %
% 'val'        - LOGICAL                                                  %
%              - Logical answer to whether the interrogated string is a   %
%                form of 'yes' or not.                                    %
% ----------------------------------------------------------------------- %
% 'word'       - STRING                                                   %
%              - Word being interrogated for a form of 'yes'.             %
%              - The following words will return true: 'Baleh', 'Da',     %
%               'Etiam', 'Go', 'Hai', 'Ja', 'K', 'Ok', 'Okay', 'Oui',     %
%               'Please', 'Si', 'Sim', 'Sure', 'Y', 'Ya', 'Ye', 'Yeah',   %
%               'Yeh', 'Yes', 'Yup'                                       %
%              - The function is not case sensitive.                      %
% ----------------------------------------------------------------------- %
%                                                                         %
% EXAMPLE                                                                 %
%                                                                         %
% Add two numbers if a governing variable is set to a form of 'yes'.      %
%                                                                         %
% >> addNums = 'Y';                                                       %
% >> if isyes(addNums)                                                    %
% disp(2 + 5);                                                            %
% end                                                                     %
%                                                                         %
% REQUIRES                                                                %
%                                                                         %
% N/A                                                                     %
%                                                                         %
% CALLED IN:                                                              %
%                                                                         %
% N/A                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% isyes
function [val] = isyes(word)

%% Define words that qualify as a form of yes.
chk = {'Baleh', 'Da', 'Etiam', 'Go', 'Hai', 'Ja', 'K', 'Ok', 'Okay',    ...
       'Oui', 'Please', 'Si', 'Sim', 'Sure', 'Y', 'Ya', 'Ye', 'Yeah',   ...
       'Yeh', 'Yes', 'Yup'}.';

%% Check if 'word' qualifies as a yes.
val = max(strcmpi(word, chk));

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
