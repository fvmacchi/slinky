function [K, C] = assembleMatrices(NODES, CONNECTIONS, k_constant, DAMPENING)

    nodesSize = size(NODES);
    conSize = size(CONNECTIONS);

    % % % Global stiffness, damper, and mass matrices
    K = zeros(nodesSize(1)*3);
    C = zeros(nodesSize(1)*3);

    % % % Iterate through list of connections. Perform local analysis
    % % % and assembly on each run through.
    for i = 1:conSize(1)
        conn = CONNECTIONS(i,:);
        local1 = NODES(conn(1),:);
        local2 = NODES(conn(2),:);
        x = local1(1) - local2(1);
        y = local1(2) - local2(2);
        z = local1(3) - local2(3);
        hyp = sqrt(x^2 + y^2 + z^2);
        x = x/hyp;
        y = y/hyp;
        z = z/hyp;
        T = [
            x^2 x*y x*z -x^2 -x*y -x*z;
            x*y y^2 y*z -x*y -y^2 -y*z;
            x*z y*z z^2 -x*z -y*z -z^2;
            -x^2 -x*y -x*z x^2 x*y x*z;
            -x*y -y^2 -y*z x*y y^2 y*z;
            -x*z -y*z -z^2 x*z y*z z^2;
        ];
        k = k_constant(i)/hyp*T;
        c = DAMPENING*T;
        K(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) = K(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) + k(1:3,1:3);
        K(conn(1)*3-2:conn(1)*3, conn(2)*3-2:conn(2)*3) = K(conn(1)*3-2:conn(1)*3, conn(2)*3-2:conn(2)*3) + k(1:3,4:6);
        K(conn(2)*3-2:conn(2)*3, conn(1)*3-2:conn(1)*3) = K(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) + k(4:6,1:3);
        K(conn(2)*3-2:conn(2)*3, conn(2)*3-2:conn(2)*3) = K(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) + k(4:6,4:6);
        C(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) = C(conn(1)*3-2:conn(1)*3, conn(1)*3-2:conn(1)*3) + c(1:3,1:3);
        C(conn(1)*3-2:conn(1)*3, conn(2)*3-2:conn(2)*3) = C(conn(1)*3-2:conn(1)*3, conn(2)*3-2:conn(2)*3) + c(1:3,4:6);
        C(conn(2)*3-2:conn(2)*3, conn(1)*3-2:conn(1)*3) = C(conn(2)*3-2:conn(2)*3, conn(1)*3-2:conn(1)*3) + c(4:6,1:3);
        C(conn(2)*3-2:conn(2)*3, conn(2)*3-2:conn(2)*3) = C(conn(2)*3-2:conn(2)*3, conn(2)*3-2:conn(2)*3) + c(4:6,4:6);
    end
end

