function output = E2MOTEndpoint(input)

%-------------------------------------------------------%
% Input:                                                %
%   auxdata                                             %
%   input.phase - a structure containing:               %
%       initialtime                                     %
%       initialstate                                    %
%       finaltime                                       %
%       finalstate                                      %
%                                                       %
% Output:                                               %
%   output - a structure containing the following:      %
%       objective                                       %
%       event                                           %
%-------------------------------------------------------%

%setting event constraints
ff = input.phase.finalstate(2);
gf = input.phase.finalstate(3);
hf = input.phase.finalstate(4);
kf = input.phase.finalstate(5);
L  = input.phase.finalstate(6);
Lm = input.phase.finalstate(9);
event1 = ff^2 + gf^2;
event2 = hf^2 + kf^2;
event3 = L - Lm;
output.eventgroup(1).event = [event1,event2,event3];

%computing the minimal time required
mf               = input.phase.finalstate(end,7);
output.objective = -mf;
