%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Combined Transfer from Insertion Orbit, Target Rendezvous, and Deorbit
% Author: Travis Hastreiter 
% Created On: 13 March, 2026
% Description: Orbit transfer using Q-Law and rendezvous using SCP in Hill 
% frame then deorbit (safely).
% Most Recent Change: 31 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Satellite Scenario
% Specify scenario
startTime = datetime(2020,5,11,12,35,38);
stopTime = startTime + days(200);
sampleTime = 60*2;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Load Satellite info
tleFile = "eccentricOrbitSatellite.tle";

satTwoBodyKeplerian = satellite(sc,tleFile, ...
    "Name","satTwoBodyKeplerian", ...
    "OrbitPropagator","two-body-keplerian");

% Specify parking orbit parameters

%% Transfer and Rendezvous from Parking Orbit to Debris
% Compute QLaw transfer from parking to debris orbit without perturbations
%QLaw_transfer

%% 

%% Optimize 

