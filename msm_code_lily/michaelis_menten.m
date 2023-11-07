% Michaelis-Menten equation. beta is a vector of Vmax and Km
function y = michaelis_menten(beta,S)

Vmax = beta(1);
Km = beta(2);
y = Vmax*S ./ (Km+S);

end