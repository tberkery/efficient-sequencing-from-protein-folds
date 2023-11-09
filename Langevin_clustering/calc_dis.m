
% Measures distance between p1 and p2
function distance = calc_dis(p1,p2)

dx = p1(1)-p2(1);
dy = p1(2)-p2(2);
dz = p1(3)-p2(3);

distance = sqrt(dx^2 + dy^2 + dz^2);

end