function[p]=probability_from_transitions(x)
    [num_rows, num_cols] = size(x);
    if (num_rows ~= num_cols)
        error("Error: Invalid transition matrix")
    end
    p = zeros(num_rows, 1);
    x_sum = sum(x, "all");
    for i = 1:num_rows
        p(i) = x(i,i)/x_sum;
    end
end