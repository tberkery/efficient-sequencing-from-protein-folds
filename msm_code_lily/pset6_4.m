% pset6_4

%% 4 a)

% Read data from file
table = readtable("onesequence_-22.79.dat");

% Get hydrophobicity table and convert to an array
hptable = table(:,1);
hp = cell2mat(table2cell(hptable));

% Get coordinates of particles table and convert to an array
coordstable = table(:,2:4);
coords = cell2mat(table2cell(coordstable));

% Implement function to find local min
[min_pos, min_E] = steepest_descent(coords,hp,10^(-6),0.01,0.5,1.01);

%% 4b)
% The criteria I used is that the potential energy difference of the new 
% and old state gets lower than 10^(-6)

%% 4c)
min_E

%% 4d)

plot3(min_pos(:,1),min_pos(:,2),min_pos(:,3),".-");
hold on
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Configuration of a 27 bead protein structure")
plot3(coords(:,1),coords(:,2),coords(:,3),"o-");


