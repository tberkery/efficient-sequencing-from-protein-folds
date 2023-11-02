function[distance]=distance(a, b)
    distance = sqrt(sum((a-b).^2));
end