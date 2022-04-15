function output = E2MOTContinuous(input)

%-------------------------------------------------------%
% Input:                                                %
%   auxdata - a structure containing the following:     %
%       time                                            %
%       state                                           %
%       control                                         %
%                                                       %
% Output:                                               %
%   output - a structure containing the following:      %
%       dynamics                                        %
%       path                                            %
%       integrand                                       %
%-------------------------------------------------------%

mu        = input.auxdata.mu;
ve        = input.auxdata.ve;
ne        = input.auxdata.ne;
nm        = input.auxdata.nm;

%calling the initial guesses for each of the following:
%  %
t           = input.phase.time(:,1);
p           = input.phase.state(:,1);
f           = input.phase.state(:,2);
g           = input.phase.state(:,3);
h           = input.phase.state(:,4);
k           = input.phase.state(:,5);
L           = input.phase.state(:,6);
m           = input.phase.state(:,7);
Le          = input.phase.state(:,8);
Lm          = input.phase.state(:,9);
i1          = input.phase.control(:,1);
i2          = input.phase.control(:,2);
i3          = input.phase.control(:,3);
T           = input.phase.control(:,4);

w           = 1 + (f.*cos(L)) + (g.*sin(L));
s           = sqrt(1 + (h.^2) + (k.^2));
deltar      = (T./m).*i1;
deltat      = (T./m).*i2;
deltan      = (T./m).*i3;

%calculating the values of the diff. equations:
pdot    = (2*p./w).*(sqrt(p/mu)).*deltat;
fdot    = sqrt(p/mu).*(deltar.*sin(L)+((w+1).*cos(L)+f).*(deltat./w)-(h.*sin(L)-k.*cos(L)).*(g.*deltan./w));
gdot    = sqrt(p/mu).*(-deltar.*cos(L)+((w+1).*sin(L)+g).*(deltat./w)+(h.*sin(L)-k.*cos(L)).*(f.*deltan./w));
hdot    = sqrt(p/mu).*((s.^2.*deltan)./(2*w)).*cos(L);     
kdot    = sqrt(p/mu).*((s.^2.*deltan)./(2*w)).*sin(L);    
Ldot    = sqrt(mu*p).*(w./p).^2 + (1./w).*sqrt(p/mu).*(h.*sin(L) - k.*cos(L)).*deltan;
mdot    = -T/ve;
Ledot	= ne*ones(length(Ldot),1);
Lmdot   = nm*ones(length(Ldot),1);

%updated guess as function is running
output(1).dynamics  = [pdot,fdot,gdot,hdot,kdot,Ldot,mdot,Ledot,Lmdot];
output(1).path      = i1.^2 + i2.^2 + i3.^2;
