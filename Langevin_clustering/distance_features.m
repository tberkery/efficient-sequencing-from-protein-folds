function[distance]=distance_features(a, b)
    distance = sqrt(sum((a-b).^2));
end