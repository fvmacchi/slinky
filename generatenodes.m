function [ NODES, CONNECTIONS, NODES_SIZE, NUMBER_OF_CONNECTIONS, NUMBER_OF_SEGMENTS ] = generatenodes( FRAMES_PER_COIL, PITCH, COIL_DIAMETER, WIRE_THICKNESS, WIRE_HEIGHT, NUMBER_OF_COILS )

theta = linspace(0, 2*pi, (FRAMES_PER_COIL+1));
NODES = [0 0 0];
COIL_RADIUS = COIL_DIAMETER/2;
INNER_COIL_RADIUS = COIL_RADIUS - WIRE_THICKNESS/2;
OUTER_COIL_RADIUS = COIL_RADIUS + WIRE_THICKNESS/2;

% % % Generate Matrix of Nodes
for i =1:NUMBER_OF_COILS
    for j = 1:FRAMES_PER_COIL
        %%%x = COIL_DIAMETER/2*cos(theta(j));
        %%%y = COIL_DIAMETER/2*sin(theta(j));
        %%%z = (i-1)*PITCH+(j-1)/FRAMES_PER_COIL*PITCH;
        
        x1 = OUTER_COIL_RADIUS*cos(theta(j));
        y1 = OUTER_COIL_RADIUS*sin(theta(j));
        z1 = ((i-1)*PITCH+(j-1)/FRAMES_PER_COIL*PITCH) + WIRE_HEIGHT/2;
        
        x2 = OUTER_COIL_RADIUS*cos(theta(j));
        y2 = OUTER_COIL_RADIUS*sin(theta(j));
        z2 = ((i-1)*PITCH+(j-1)/FRAMES_PER_COIL*PITCH) - WIRE_HEIGHT/2;
        
        x3 = INNER_COIL_RADIUS*cos(theta(j));
        y3 = INNER_COIL_RADIUS*sin(theta(j));
        z3 = ((i-1)*PITCH+(j-1)/FRAMES_PER_COIL*PITCH) + WIRE_HEIGHT/2;
        
        x4 = INNER_COIL_RADIUS*cos(theta(j));
        y4 = INNER_COIL_RADIUS*sin(theta(j));
        z4 = ((i-1)*PITCH+(j-1)/FRAMES_PER_COIL*PITCH) - WIRE_HEIGHT/2;
        
        x5 = COIL_RADIUS*cos((theta(j)+theta(j+1))/2);
        y5 = COIL_RADIUS*sin((theta(j)+theta(j+1))/2);
        z5 = ((i-1)*PITCH+(j-0.5)/FRAMES_PER_COIL*PITCH);
        
        
        NODES = [NODES;
                    x1 y1 z1;
                    x2 y2 z2;
                    x3 y3 z3;
                    x4 y4 z4;
                    x5 y5 z5;];
    end
end

% % % Remove dummy initial node [0 0 0]
NODES = NODES(2:size(NODES,1),:);
NODES_SIZE = size(NODES,1);

% % % Have slinky coil downwards instead
NODES(:,3) = -1*NODES(:,3);

NUMBER_OF_SEGMENTS = NUMBER_OF_COILS*FRAMES_PER_COIL;

% % % Generate matrix of connections
NUMBER_OF_CONNECTIONS = NUMBER_OF_SEGMENTS*16 - 8;
CONNECTIONS = zeros(NUMBER_OF_CONNECTIONS,2);

for i = 1:(NUMBER_OF_SEGMENTS)
    frame_a_node_1 = (i-1)*5 + 1;
    frame_a_node_2 = (i-1)*5 + 2;
    frame_a_node_3 = (i-1)*5 + 3;
    frame_a_node_4 = (i-1)*5 + 4;
    frame_a_node_5 = (i-1)*5 + 5;
    
    frame_b_node_1 = (i-1)*5 + 6;
    frame_b_node_2 = (i-1)*5 + 7;
    frame_b_node_3 = (i-1)*5 + 8;
    frame_b_node_4 = (i-1)*5 + 9;
    
    CONNECTIONS((i-1)*16+1,:) = [frame_a_node_1 frame_a_node_2];
    CONNECTIONS((i-1)*16+2,:) = [frame_a_node_2 frame_a_node_3];
    CONNECTIONS((i-1)*16+3,:) = [frame_a_node_3 frame_a_node_4];
    CONNECTIONS((i-1)*16+4,:) = [frame_a_node_4 frame_a_node_1];
    
    CONNECTIONS((i-1)*16+5,:) = [frame_a_node_1 frame_a_node_5];
    CONNECTIONS((i-1)*16+6,:) = [frame_a_node_2 frame_a_node_5];
    CONNECTIONS((i-1)*16+7,:) = [frame_a_node_3 frame_a_node_5];
    CONNECTIONS((i-1)*16+8,:) = [frame_a_node_4 frame_a_node_5];
    
    if i ~= NUMBER_OF_SEGMENTS
        CONNECTIONS((i-1)*16+9,:) = [frame_a_node_5 frame_b_node_1];
        CONNECTIONS((i-1)*16+10,:) = [frame_a_node_5 frame_b_node_2];
        CONNECTIONS((i-1)*16+11,:) = [frame_a_node_5 frame_b_node_3];
        CONNECTIONS((i-1)*16+12,:) = [frame_a_node_5 frame_b_node_4];

        CONNECTIONS((i-1)*16+13,:) = [frame_a_node_1 frame_b_node_1];
        CONNECTIONS((i-1)*16+14,:) = [frame_a_node_2 frame_b_node_2];
        CONNECTIONS((i-1)*16+15,:) = [frame_a_node_3 frame_b_node_3];
        CONNECTIONS((i-1)*16+16,:) = [frame_a_node_4 frame_b_node_4];
    end
end

end

