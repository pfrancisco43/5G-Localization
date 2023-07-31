% Author: Paulo Francisco
% 
% This script provides the necessary parameters to test the geometry-based 
% localization method.
% 
% If you change the positioning of the BS, MS, or SC, please also modify 
% the search grid in the parameterEstimation.m file.

clc;clear;close all
%% INPUTS
los=1;              % if los in this simulation (0-no, 1-yes)
b=[-8, 0, 5]';      % BS
m=[7, 10, 1]';      % MS 
s1=[-10, 4, 3]';    % SC_1
s2=[-10, 8, 4]';    % SC_2
s3=[10, 8, 3]';     % SC_3
s4=[10, 4, 4]';     % SC_4
s=[s1 s2 s3 s4];    % set of scatterers
Tol=0.5;            % tolerance value for considering the existence of a LoS path
c=300;              % Propagation speed in us
Nt=32;              % number of TX antennas
Nr=32;              % number of RX antennas
N=10;               % number of subcarriers
B=100;              % total BW in MHz
Ns=20;              % number of sybols sent
L_az=14;            % grid resolution for azimuth
L_el=14;            % grid resolution for elevation
fin=1;              % Fining step in the DCS-SOMP

%rng(1)
%% Get True parameters
[ToAs, AoDs, AoAs, L]=getTrueParameters(b,m,s,los,c);

%% Stage 1 - Channel modeling 
[y,x,H]=channelModeling (ToAs, AoDs, AoAs, Nt, Nr, N, B, Ns);

%% Stage 2 - Parameter Estimation
[ToAs_E, AoDs_E, AoAs_E]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, H, fin);

%% Verify if LoS exists
[haveLos]=verifyLos(AoDs_E, AoAs_E, Tol);

%% Stage 3 - geometry-based algorithms
if haveLos    
    [err,pos]=LoSAlgorithm(b, ToAs_E(1), AoDs_E(1,:), m, c);
    fprintf('LoS exists: [YES]\nError: %.2f\nEstimated MS postion: (%.1f %.1f %.1f)\n',err,pos);
    showLoSScenario(b,m,pos);
else
    [err, pos_m, pos_s, u, v, k]=NLoSAlgorithm(b, ToAs_E, AoDs_E, AoAs_E, m, c);
    fprintf('LoS exists: [NO]\nError: %.2f\nEstimated MS postion: (%.1f %.1f %.1f)\n',err,pos_m);
    showNLoSScenario(b,m,s,pos_m,pos_s,u,v,k)
end
