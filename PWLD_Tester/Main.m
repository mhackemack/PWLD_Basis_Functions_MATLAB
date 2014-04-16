%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D/3D PWLD Run Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% -------------------
clc;clear all;close all;
format long e
% Populate path with additional folders
% -------------------------------------
addpath(genpath('..'));   % add all folders/subfolders
global glob
glob = setGlobals();
% Load data and execute program
% -----------------------------
[data, geometry] = load_user_input();
[sol, geometry, DoF] = execute_problem(data, geometry);

