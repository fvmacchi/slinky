function [I] = LUDecomp(a, b)
    array_order = size(a,1);
    L = zeros(array_order(1));
    U = zeros(array_order);
    
    for i = 1:array_order
        U(i,i) = 1; 
    end
    % Step 1
    for i = 1:array_order
        L(i, 1) = a(i, 1);
    end
    % Step 2
    for j = 2:array_order
       U(1,j) = a(1,j)/L(1,1);
    end
    % Step 3
    for j = 2:array_order-1
       % Part 1
       for i = j:array_order
           L(i,j) = a(i,j);
           for k=1:j-1
               L(i,j) = L(i,j) - (L(i,k)*U(k,j));
           end
       end
       % part 2
       for k = j+1:array_order
           U(j,k) = a(j,k);
           for i = 1:j-1
               U(j,k) = U(j,k) - L(j,i) * U(i,k);
           end
           U(j,k) = U(j,k)/L(j,j);
       end
    end

    % Step 4
    L(array_order,array_order) = a(array_order,array_order);
    for k = 1:array_order-1
        L(array_order,array_order) = L(array_order,array_order) - L(array_order,k) * U(k,array_order);
    end
    I = zeros(array_order);
    
    D = zeros(array_order,1);
    X = zeros(array_order,1);
    for i = 1:array_order
        D(i) = b(i);
        for k = 1:i-1
            D(i) = D(i) - L(i,k)*D(k);
        end
        D(i) = D(i) / L(i,i);
    end
    for i=array_order:-1:1
        X(i) = D(i);
        for k = array_order:-1:i+1
            X(i) = X(i) - U(i,k)*X(k);
        end
    end
    I = X;
end