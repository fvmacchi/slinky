function [ANSWER] = gaussseidel (A, B)
    array_order = size(A,1);
    ANSWER = zeros(array_order,1);
    
    NewValues = zeros(array_order,1);
    OldValues = zeros(array_order,1);
    % Set OldValues to have non zero number, so while loop will work
    OldValues(1) = 10000;
    %abs(sum(NewValues) - sum(OldValues))

    while(abs(sum(NewValues) - sum(OldValues)) > 0.0005*array_order)
    %for k = 1:20
        OldValues(:,1) = NewValues(:,1);
        for i=1:array_order
            NewValues(i,1) = B(i);
            for j = 1:i-1
                NewValues(i,1) = NewValues(i,1) - A(i,j) * NewValues(j);
            end
            for j = i+1:array_order
                NewValues(i,1) = NewValues(i,1) - A(i,j) * NewValues(j);
            end 
            NewValues(i,1) = NewValues(i,1) / A(i,i);
        end
    end
    ANSWER(:,1) = NewValues(:,1)
end
