function distance = calc_dis(x1,y1,z1,x2,y2,z2)

dx = x1-x2;
dy = y1-y2;
dz = z1-z2;

distance = sqrt(dx^2 + dy^2 + dz^2);

end