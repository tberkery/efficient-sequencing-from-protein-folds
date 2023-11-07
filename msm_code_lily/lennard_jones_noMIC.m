% function that calculates the lennard jones potential between 2 particles
%dx, dy, and dz inputs are MIC transformed
% Requires MIC.m file 
function potential = lennard_jones_noMIC(x1,y1,z1,x2,y2,z2,epsilon,sigma)

%calculate distance rij = |ri-rj|
xij = x1-x2;
yij = y1-y2;
zij = z1-z2;
rij = sqrt(xij^2 + yij^2 + zij^2);

%calculate potential step by step using the lennard jones equation
inner_bracket = sigma/rij;
outer_bracket = inner_bracket.^12 - inner_bracket.^6;
potential = 4 * epsilon * outer_bracket;
return