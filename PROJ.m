%{
GROUP 13
PROJECT 2B

NOTE: all units are in their base SI form
%}

% % % initialization
clear all;
close all;
clc;

% % % INPUTS FOR DEFINING INITIAL STATE OF SLINKY
% % % cross sectional slices in one turn
FRAMES_PER_COIL = 16;
% % % distance between centrelines of adjacent coils when compressed
PITCH = 0.00065;
COIL_DIAMETER = (0.06349 + 0.06849)/2;
% % % wire guage
WIRE_THICKNESS = 0.0025;
WIRE_HEIGHT = PITCH;
NUMBER_OF_COILS = 2;
TOTAL_SLINKY_MASS = 0.198*(NUMBER_OF_COILS/80);

GRAVITY = -9.81;
ATTACHED_MASS = 1;

% % % Find K value
E = 16.7e6;
A = WIRE_THICKNESS*WIRE_HEIGHT;
L = pi*COIL_DIAMETER/FRAMES_PER_COIL;
K0 = E*A/L;
theta = acos(L/sqrt(WIRE_THICKNESS^2 + WIRE_HEIGHT^2 + L^2));
k_constant = K0/(4+2*cos(theta));
K_constant = @(i) k_constant;

% % % Constants for implicit dynamics
BETA = 0.25;
GAMMA = 0.50;
DELTA_T = 0.1;
TIME_LIMIT = 0.5;

% % % Create node an connection matrices
[NODES, CONNECTIONS, NODES_SIZE, NUMBER_OF_CONNECTIONS, NUMBER_OF_SEGMENTS] = generatenodes(FRAMES_PER_COIL, PITCH, COIL_DIAMETER, WIRE_THICKNESS, WIRE_HEIGHT, NUMBER_OF_COILS);

% % % Genereate global mass matrix, GM
GM = zeros(NODES_SIZE*3);
MASS_PER_SEGMENT = TOTAL_SLINKY_MASS/NUMBER_OF_SEGMENTS;
sample_segment_lengths = zeros(20,1);
total_segment_length = 0;
for i = 1:20
    node1 = NODES(CONNECTIONS(i,1),:);
    node2 = NODES(CONNECTIONS(i,2),:);
    sample_segment_lengths(i) = sqrt((node1(1)-node2(1))^2 + (node1(2)-node2(2))^2 + (node1(3)-node2(3))^2);
    total_segment_length = total_segment_length + sample_segment_lengths(i);
end
mass_per_memeber_length = MASS_PER_SEGMENT/total_segment_length;
for i = 1:NUMBER_OF_CONNECTIONS
    % % % Formulate local mass matrix
    node1 = NODES(CONNECTIONS(i,1),:);
    node2 = NODES(CONNECTIONS(i,2),:);
    connection_length = sqrt((node1(1)-node2(1))^2 + (node1(2)-node2(2))^2 + (node1(3)-node2(3))^2);
    local_m = connection_length*mass_per_memeber_length/2*eye(6);
    % % % Assemble global mass matrix
    %%%display(GM(3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1), 3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1)));
    GM(3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1), 3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1)) = GM(3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1), 3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1)) + local_m(1:3,1:3);
    GM(3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1), 3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2)) = GM(3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1), 3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2)) + local_m(1:3,4:6);
    GM(3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2), 3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1)) = GM(3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2), 3*CONNECTIONS(i,1)-2:3*CONNECTIONS(i,1)) + local_m(4:6,1:3);
    GM(3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2), 3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2)) = GM(3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2), 3*CONNECTIONS(i,2)-2:3*CONNECTIONS(i,2)) + local_m(4:6,4:6);
end
GM(NODES_SIZE*3-2,NODES_SIZE*3-2) = GM(NODES_SIZE*3-2,NODES_SIZE*3-2) + ATTACHED_MASS;
GM(NODES_SIZE*3-1,NODES_SIZE*3-1) = GM(NODES_SIZE*3-1,NODES_SIZE*3-1) + ATTACHED_MASS;
GM(NODES_SIZE*3,NODES_SIZE*3) = GM(NODES_SIZE*3,NODES_SIZE*3) + ATTACHED_MASS;

% % % Used for visualizing. Very expensive to run
if 1
    figure(1);
    scatter3(NODES(:,1),NODES(:,2),NODES(:,3));
    hold on;
    for i=1:NUMBER_OF_CONNECTIONS
            tmp = [NODES(CONNECTIONS(i,1),:) ; NODES(CONNECTIONS(i,2),:)];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3), '-g');
            hold on;
    end
end

% % % Initialize vectors for displacement, velocity, acceleration and force
time = 0;
U0 = zeros((3*NODES_SIZE), 1);
V0 = zeros((3*NODES_SIZE), 1);
A0 = zeros((3*NODES_SIZE), 1);
F0 = zeros((3*NODES_SIZE), 1);
U1 = zeros((3*NODES_SIZE), 1);
V1 = zeros((3*NODES_SIZE), 1);
A1 = zeros((3*NODES_SIZE), 1);
F = zeros((3*NODES_SIZE), 1);

for i=5:NODES_SIZE
   A0(i*3) = GRAVITY;
end

F = zeros(NODES_SIZE*3, 1);
for i = 1:NODES_SIZE
    F(i*3) = GM(i*3,i*3)*GRAVITY;
end

% % % Iterate through timesteps
while time < TIME_LIMIT
    % % % FRANCESCO: plug in your bits/parts/pieces here (whose whichever sounds least dirty)
    [GK, GC] = assembleMatrices(NODES, CONNECTIONS, K_constant, 0);

    
    helper_A = 2/(BETA*DELTA_T^2)*GM + GK + (2*GAMMA)/(BETA*DELTA_T)*GC;
    helper_B = 2/(BETA*DELTA_T^2)*GM + (2*GAMMA)/(BETA*DELTA_T)*GC;
    helper_C = 2/(BETA*DELTA_T)*GM + (2*GAMMA/BETA - 1)*GC;
    helper_D = (1 - BETA)/BETA*GM + DELTA_T*GC*((GAMMA - 1)+(1 - BETA)/BETA*GAMMA);   

    % % % Use gauss siedel to solve for U1 and F1
    [U1, INTERNAL_F] = solver(helper_A, U1, F + helper_B*U0 + helper_C*V0 + helper_D*A0, [1:12], [13:NODES_SIZE*3]);
    
    % % % determine the next snapshot of the acceleration vector
    A1 = 2/(BETA*DELTA_T)*(U1-U0)/DELTA_T - 2/(BETA*DELTA_T)*V0 - (1 - BETA)/BETA*A0;

    % % % determine the next snapshot of the velocity vector
    V1 = DELTA_T*((1-GAMMA)*A0 - GAMMA*A1) + V0;

    % % % Iterate to next time
    time = time + DELTA_T; 
    U0 = U1;
    V0 = V1;
    A0 = A1;
end

% % % Post visualization
if 1
    NEW_NODES = NODES;
    for i=1:NODES_SIZE
       NEW_NODES(i,1) = NEW_NODES(i,1) + U0(i*3-2);
       NEW_NODES(i,2) = NEW_NODES(i,2) + U0(i*3-1);
       NEW_NODES(i,3) = NEW_NODES(i,3) + U0(i*3);
    end
  
    figure(2);
    scatter3(NEW_NODES(:,1),NEW_NODES(:,2),NEW_NODES(:,3));
    hold on;
    for i=1:NUMBER_OF_CONNECTIONS
            tmp = [NEW_NODES(CONNECTIONS(i,1),:) ; NEW_NODES(CONNECTIONS(i,2),:)];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3), '-g');
            hold on;
    end
end
