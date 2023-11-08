%% 1a)
load rate_vs_conc.dat
conc = rate_vs_conc(:,1);
rate = rate_vs_conc(:,2);

%% 1b) and c)
[params,~,~,COVB] = nlinfit(conc,rate,@michaelis_menten,[2,1])
params

plot(conc,rate, ".")
hold on
curve_conc = linspace(0,5,1000);
curve_rate = michaelis_menten(params,curve_conc);
plot(curve_conc,curve_rate);
xlabel("[S]")
ylabel("Reaction rate")
legend("Data points", "Fitted curve")
title("Michaelis-Menten")
hold off

%% 1d)
variance_Vmax = COVB(1,1)
variance_Km = COVB(2,2)
