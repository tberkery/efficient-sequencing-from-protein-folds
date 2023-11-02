function[x_norm_by_row]=norm_rows(x)
    [num_rows, num_cols] = size(x);
    x_norm_by_row = zeros(num_rows, num_cols);
    for i = 1:num_rows
        row_sum = sum( x(i,:) );
        x_norm_by_row(i,:) = x(i,:)/row_sum;
    end
end