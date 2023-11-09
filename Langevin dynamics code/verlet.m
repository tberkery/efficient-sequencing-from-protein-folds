% Verlet used to update position of particles for MD
function x_new = verlet(x_old, v_old, Fx_old, m, dt)
x_new = x_old + (v_old.*dt) + (dt.^2).*(Fx_old/(2.*m));
end