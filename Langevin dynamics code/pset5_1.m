% Pset 5 question 1 Langevin dynamics
clc; clear;
rng(5)

R_debug = 0; % Global variable for debugging. See random_force for usage
init_pos = [0,0,1; 0,0,2; 0,0,3; 0,0,4; 0,0,5];

% running LD with T=2, dt=0.003, steps = 10,000 (need to change)
[times, potentials, kinetics, temperatures, eq_pos] = LD(init_pos, 2, 0.003, 10000);

int_energies = kinetics + potentials;
% figure(1);
hold on
plot(times, potentials);
plot(times, kinetics);
%plot(times, temperatures)
legend("PE", "KE")
xlabel("Time")
ylabel("Energy")
hold off

ave_PE = sum(potentials(9001:10000))/length(potentials(9001:10000))
ave_KE = sum(kinetics(9001:10000))/length(kinetics(9001:10000)) 
ave_E = sum(int_energies(9001:10000))/length(int_energies(9001:10000))
%plot3(eq_pos(:,1), eq_pos(:,2), eq_pos(:,3))
