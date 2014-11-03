clear all
% 'Biomechanics and Motor Control of Human Movement' by David Winter (4th edition)
% Sample kinematic data set
% saved as TEST.m
 
% kinetic data from pg 191
% Reaction forces (RD from frame 7)  and Moment forces (M from frame 7)
RD7 = [-134.33 -768.38 -10.26]; % units are Newtons in GRS
MD7 = [16.13 -0.49 -101.32]  ;   % units are Newtons*meters in GRS
% RP7 = [-96.41 -755.17 -13.15]; % proximal reaction force used in the g to a axes system 
% kinetic data 'table 7.4' Leg angular displacements
theta7 = [-0.14815 -0.02705 -0.51094]; %units are radians
theta6 = [-0.14781 -0.02501 -0.46053];
theta8 = [-0.15161 -0.02809 -0.56749];
% data from table 7.4, angular velocity. units is rad/s
thetadot7 = [-0.1139 -0.0923 -3.209];
thetadot6 = [-0.0911 -0.1048 -2.909];
thetadot8 = [-0.1109 0.2175 -3.472];
 
ax = 7.029;
ay = 1.45;
az = -0.348;
ix = 0.0138; % units are kg*m^2
iy = 0.0024; % units are kg*m^2
iz = 0.0138; % units are kg*m^2
ld = 0.1386; % in meters
lp = 0.1815; % in meters
g = 9.81; % gravity in m/s^2
M = 3.22; % mass in kg
pi = 3.14;% constant pi
 
% Step 1, RP is proxial, RD is distanl
RP7(1)= RD7(1) + (M*ax);
RP7(2)= RD7(2) + (M*g)+ (M*ay);
RP7(3)= RD7(3) + (M*az);
 
%matrix equation,  commas indicate new line, while spaces indicate new term
matrix = [ cos(theta7(2))*cos(theta7(3)) sin(theta7(3))*cos(theta7(1))+sin(theta7(1))*sin(theta7(2))*cos(theta7(3)) sin(theta7(1))*sin(theta7(3))-cos(theta7(1))*sin(theta7(2))*cos(theta7(3)); -(cos(theta7(2))*sin(theta7(3))) (cos(theta7(1))*cos(theta7(3))-sin(theta7(1))*sin(theta7(2))*sin(theta7(3))) (sin(theta7(1))*cos(theta7(3))+cos(theta7(1))*sin(theta7(2))*sin(theta7(3))); sin(theta7(2)) -(sin(theta7(1))*cos(theta7(2))) (cos(theta7(1))*(cos(theta7(2))))];
% GTA = [0.8718 -0.4800 0.0954, 0.4887 0.8645 -0.1156, -0.027 0.1475 0.9886];
% need to transform the reaction forces from the global to the anatomical
% axes system
GTA = [0.8718 -0.4800 0.0954, 0.4887 0.8645 -0.1156, -0.027 0.1475 0.9886];
% GTA stands for the global to anatomical matrix
Rd7 = (matrix)*(RD7)';
Rp7 = (matrix)*(RP7)';
 
%step 3, transform ankle moments from g to a axes system
Md7 = (matrix)*(MD7)';
 
%calculate angular velocities and acceleration for usage in Euler's kinetic
%equation (7.9)
% wx is the first term, wy is the second, and wz is the third term.
w = [(cos(theta7(2))*cos(theta7(3))) (sin(theta7(3))) (0);  -cos(theta7(2))*sin(theta7(3)) (cos(theta7(3))) (0); sin(theta7(2)) (0) (1)]
w7 = w*thetadot7' 
w6 = w*thetadot6'
w8 = w*thetadot8'
 
% to calculate the angular accelerations of the segment where delta t is
% the sampling period....to get the difference of 7, 8 and 6 must be
% utilized
dt = 0.03333; % this is delta of time
alp7  = [w8 - w6] / 2*dt % answer is in r/s^2
 
% solve Eulers equations (7.9) for the proximal moments 
 
Mxp =(ix)*(alp7(1)) + (iz-iy)*(w(2)*w(3)) -Rd7(3)*ld -Rp7(3)*lp + Md7(1)  % units are Newton*Meters
Myp = (iy)*(alp7(2)) + (ix-iz)*(w(1)*w(3)) + Md7(2)
Mzp = (iz)*(alp7(3)) + (iy-ix)*(w(1)*w(2))+ Rd7(1)*ld + Rp7(1)*lp + Md7(3)
