%{
GROUP 13
PROJECT 2B

NOTE: all units are in their base SI form
%}

% % % initialization
clear all;
close all;
clc;

E = 16.7e6;

% % % INPUTS FOR DEFINING INITIAL STATE OF SLINKY
% % % cross sectional slices in one turn
FRAMES_PER_COIL = 16;
% % % distance between centrelines of adjacent coils when compressed
PITCH = 0.00065;
% % % distance between centrelines of adjacent coils at yield
YIELD_PITCH = 0.084;
COIL_DIAMETER = (0.06349 + 0.06849)/2;
% % % wire guage
WIRE_THICKNESS = 0.0025;
WIRE_HEIGHT = 0.00065;
NUMBER_OF_COILS = 2;

[NODES, CONNECTIONS, NODES_SIZE, NUMBER_OF_CONNECTIONS] = generatenodes(FRAMES_PER_COIL, PITCH, COIL_DIAMETER, WIRE_THICKNESS, WIRE_HEIGHT, NUMBER_OF_COILS);
[YIELD_NODES, YIELD_CONNECTIONS, NODES_SIZE, NUMBER_OF_CONNECTIONS] = generatenodes(FRAMES_PER_COIL, YIELD_PITCH, COIL_DIAMETER, WIRE_THICKNESS, WIRE_HEIGHT, NUMBER_OF_COILS);



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
% % % Used for visualizing. Very expensive to run
if 1
    figure(2);
    scatter3(YIELD_NODES(:,1),YIELD_NODES(:,2),YIELD_NODES(:,3));
    hold on;
    for i=1:NUMBER_OF_CONNECTIONS
            tmp = [YIELD_NODES(CONNECTIONS(i,1),:) ; YIELD_NODES(CONNECTIONS(i,2),:)];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3), '-g');
            hold on;
    end
end

STRESSES = zeros(NUMBER_OF_CONNECTIONS,1);
for i = 1:NUMBER_OF_CONNECTIONS
    connection = CONNECTIONS(i,:);
    a = NODES(connection(1),:);
    b = NODES(connection(2),:);
    x = b(1) - a(1);
    y = b(2) - a(2);
    z = b(3) - a(3);
    l0 = sqrt(x^2 + y^2 + z^2);
    connection = YIELD_CONNECTIONS(i,:);
    a = YIELD_NODES(connection(1),:);
    b = YIELD_NODES(connection(2),:);
    x = b(1) - a(1);
    y = b(2) - a(2);
    z = b(3) - a(3);
    l = sqrt(x^2 + y^2 + z^2);
    strain = (l0 - l)/l0;
    STRESSES(i) = strain*E;
end
display(abs(min(STRESSES)));
    