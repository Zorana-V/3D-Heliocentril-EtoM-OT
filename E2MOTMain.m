clc; clear;

%---------------------------------------------------------------%
% Helicenric, 3D Orbit Transfer from Earth to Mars (Part 1)     %
%---------------------------------------------------------------%
% Objective: Find a 3D solution for the orbit transfer from     %
%   Earth and Mars during the heliocentic phase using modified  %
%   equinoctial elements.                                       %
%                                                               %
% Constraints:                                                  %
%   Diff Eqns of Motion:                                        %
%       pdot = sqrt(p/mu)(2p/w)*delt                            %
%       fdot = sqrt(p/mu)[(delr)sinL + ((w+1)cosL + f)(delt/w)  %
%                   - ((h)sinL - (k)cosL)(g*deln)/w]            %
%       gdot = sqrt(p/mu)[(-delr)cosL + (w+1)sinL + g)(delt/w)  %
%                   + ((h)sinL - (k)cosL)(g*deln)/w]            %
%       hdot = sqrt(p/mu)(s^2)(deln/2w)cosL                     %
%       kdot = sqrt(p/mu)(s^2)(deln/2w)sinL                     %
%       Ldot = sqrt(p*mu)(w/p)^2                                %
%                   + sqrt(p/mu)((h)sinL - (k)cosL)(deln/w)     %
%       mdot    = -T/ve                                         %
%       Ledot	= ne                                            %
%       Lmdot	= nm                                            %
%                                                               %
%   Bounds:                                                     %
%       (p(0),p(tf)) = (p0,pf)                                  %
%       (f(0),f(tf)) = (-1,1)                                   %
%       (g(0),g(tf)) = (-1,1)                                   %
%       (h(0),h(tf)) = (-1,1)                                   %
%       (k(0),k(tf)) = (-1,1)                                   %
%       (L(0),L(tf)) = (L0,Lf)                                  %
%       (m(0),m(tf)) = (m0,free)                                %
%       (Le(0),Le(tf) = (L0,free)                               %
%       (Lm(0),Lm(tf) = (free,Lf)                               %
% --------------------------------------------------------------%

%---------------------------------------------------------------%
% Defining Initial Conditions, Universal Constants, & Boundary  %
% Conditions                                                    %
%---------------------------------------------------------------%

mu              = 1;                %gravitational parameter of the Sun [km^3 * s^-1]
T               = 0.1405;           %maximum thrust
ve              = 1.87583;          %escape velocity earth


% BOUNDS ON VARIABLES
t0              = 0;                %initial time
r0              = 1;                %initial radius (Earth Radius) [AU]
a0              = r0;               %classical orbital elements earth
e0              = 0;
bOmega0         = 0;
inc0            = 0;
lomega0         = 0;
nu0             = 0;
m0              = 1;                %initial mass
ne              = 1/sqrt(a0^3);     %mean motion of earth

rf              = 1.5;              %final radius (Mars Radius) [AU]
af              = rf;               %classical orbital elements mars
ef              = 0;
bOmegaf         = 0;
incf            = 0;
lomegaf         = 0;
nuf             = 0;
nm              = 1/sqrt(af^3);     %mean motion of mars

a = (a0 + af)/2;                    %semi-major axis of orbit transfer
htau = pi*sqrt((a^3)/mu);           %half of the period of orbit transfer

coe0            = [a0; e0; bOmega0; inc0; lomega0; nu0];
coef            = [af; ef; bOmegaf; incf; lomegaf; nuf];

mee0            = coe2mee(coe0);
meef            = coe2mee(coef);

p0 = mee0(1);
f0 = mee0(2);
g0 = mee0(3);
h0 = mee0(4);
k0 = mee0(5);
L0 = mee0(6);
Le0= L0;

pf = meef(1);
ff = meef(2);
gf = meef(3);
hf = meef(4);
kf = meef(5);
Lf = meef(6);

tfmin       = t0;                   %lower tolerance on final time
tfmax       = 10;                   %upper tolerance on final time
rfmin       = r0;                   %lower tolerance on final radius
rfmax       = 10;                   %upper tolerance on final radius
pmin        = p0;                   %lower tolerance of p
pmax        = 10*p0;                %upper tolerance of p
fmin        = -1;                   %lower tolerance of f
fmax        = +1;                   %upper tolerance of f
gmin        = -1;                   %lower tolerance of g
gmax        = +1;                   %upper tolerance of g
hmin        = -10;                  %lower tolerance of h
hmax        = +10;                  %upper tolerance of h
kmin        = -10;                  %lower tolerance of k
kmax        = +10;                  %upper tolerance of k
Lmin        = 0;                    %lower tolerance of true longitude
Lmax        = 6*pi;                 %upper tolerance of true longitude
mmin        = 0.1;                  %lower tolerance on mass
mmax        = m0;                   %upper tolerance on mass
i1min       = -10;                  %lower tolerance on deltar unit vector
i1max       = +10;                  %upper tolerance on deltar unit vector
i2min       = -10;                  %lower tolerance on deltat unit vector
i2max       = +10;                  %upper tolerance on deltat unit vector
i3min       = -10;                  %lower tolerance on deltan unit vector
i3max       = +10;                  %upper tolerance on deltan unit vector
Tmin        = 0;                    %lower tolerance on thrust
Tmax        = T;                    %upper tolerance on thrust
Lemin       = 0;                    %lower tolerance of true longitude of earth
Lemax       = 6*pi;                 %upper tolerance of true longitude of earth
Lmmin       = 0;                    %lower tolerance of true longitude of mars
Lmmax       = 6*pi;                 %upper tolerance of true longitude of mars

% AUXDATA = Auxiliary Data [Structure]
auxdata.mu      = mu;
auxdata.ve      = ve;
auxdata.ne      = ne;
auxdata.nm      = nm;

% BOUNDS [STRUCTURE]
bounds.phase.initialtime.lower      = t0;
bounds.phase.initialtime.upper      = t0;
bounds.phase.finaltime.lower        = tfmin;
bounds.phase.finaltime.upper        = tfmax;
bounds.phase.initialstate.lower     = [p0,f0,g0,h0,k0,L0,m0,Le0,Lmmin];
bounds.phase.initialstate.upper     = [p0,f0,g0,h0,k0,L0,m0,Le0,Lmmax];
bounds.phase.state.lower            = [pmin,fmin,gmin,hmin,kmin,Lmin,mmin,Lemin,Lmmin];
bounds.phase.state.upper            = [pmax,fmax,gmax,hmax,kmax,Lmax,mmax,Lemax,Lmmax];
bounds.phase.finalstate.lower       = [pf,fmin,gmin,hmin,kmin,Lmin,mmin,Lemin,Lmmin];
bounds.phase.finalstate.upper       = [pf,fmax,gmax,hmax,kmax,Lmax,mmax,Lemax,Lmmax];
bounds.phase.control.lower          = [i1min,i2min,i3min,Tmax];
bounds.phase.control.upper          = [i1max,i2max,i3max,Tmax];

bounds.phase.path.lower             = 1;
bounds.phase.path.upper             = 1;

bounds.eventgroup(1).lower          = [ef^2, tan(incf/2)^2, 0];
bounds.eventgroup(1).upper          = [ef^2, tan(incf/2)^2, 0];


% GUESS [STRUCTURE]

tGuess  = [t0;tfmax];
pGuess  = [p0;pf];
fGuess  = [f0;f0];
gGuess  = [g0;g0];
hGuess  = [h0;h0];
kGuess  = [k0;k0];
LGuess  = [Lmin;Lmin];
mGuess  = [mmax;mmin];
LeGuess = [Le0;Lemax];
LmGuess = [Lmmin;Lmmax];
i1Guess = [-1; 1];
i2Guess = [-1; 1];
i3Guess = [0; 0];
TGuess = [0; 0];

guess.phase.time    = tGuess;
guess.phase.state   = [pGuess,fGuess,gGuess,hGuess,kGuess,LGuess,mGuess,LeGuess,LmGuess];
guess.phase.control = [i1Guess,i2Guess,i3Guess,TGuess];

% MESH [STRUCTURE]
numIntervals         = 10;
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-7;
mesh.maxiterations   = 5;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
mesh.phase.colpoints = 3*ones(1,numIntervals);
mesh.phase.fraction  = ones(1,numIntervals)/numIntervals;

% SETUP [STRUCTURE]
setup.name                           = 'Earth2Mars-Orbit-Transfer-Problem';
setup.functions.continuous           = @E2MOTContinuous;
setup.functions.endpoint             = @E2MOTEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-6;
setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'none';
setup.method                         = 'RPM-Differentiation';

% CALL GPOPS-II
output = gpops2(setup); %function to use gpops2 w Setup Structure

%finding vectors that represent the orbit
solution    = output.result.solution;
t           = solution.phase.time(:,1);
p           = solution.phase.state(:,1);
f           = solution.phase.state(:,2);
g           = solution.phase.state(:,3);
h           = solution.phase.state(:,4);
k           = solution.phase.state(:,5);
L           = solution.phase.state(:,6);
m           = solution.phase.state(:,7);
Le          = solution.phase.state(:,8);
Lm          = solution.phase.state(:,9);
i_1         = solution.phase.control(:,1);
i_2         = solution.phase.control(:,2);
i_3         = solution.phase.control(:,3);
T           = solution.phase.control(:,4);

Rv = zeros(length(p),3);            %intializing a position vector for orbit transfer

%finding position of orbit transfer
for ii = 1:length(p)
    mee = [p(ii); f(ii); g(ii); h(ii); k(ii); L(ii)];
    [Rv(ii,:),~] = mee2rv(mee);
end

%finding posiiton of earth and mars
[rv1,~] = mee2rv([p0 f0 g0 h0 k0 Le(1)]);
[rv2,~] = mee2rv([pf ff gf hf kf Lm(1)]);

Orbit1 = rv1;
Orbit2 = rv2;

for ii = 2:length(Lm)
    [rv1,~] = mee2rv([p0 f0 g0 h0 k0 Le(ii)]);
    [rv2,~] = mee2rv([pf ff gf hf kf Lm(ii)]);
    Orbit1 = [Orbit1; rv1];
    Orbit2 = [Orbit2; rv2];
end

%plotting orbit transfer
plot3(Rv(:,1),Rv(:,2),Rv(:,3))
hold on
plot3(Orbit1(:,1),Orbit1(:,2),Orbit1(:,3))
hold on
plot3(Orbit2(:,1),Orbit2(:,2),Orbit2(:,3))
legend('Transfer Orbit','Earth Orbit','Mars Orbit')

