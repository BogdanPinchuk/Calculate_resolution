function dovzin = dovzdotoch(h, Rk, alpha)
% Функція розрахунку відстані від КА до точки пікселя на Землі

% dovzin - відстань від КА до кожної точки яка описує РЗ пікселя
% h - скорегована висота супутника 
% Rk - радіус кривизни Землі
% alpha - кут відхилення від надиру

beta = asin(((h + Rk) ./ Rk) .* sin(alpha)) - alpha;

dovzin = (h + Rk .* (1 - cos(beta))) ./ cos(alpha);

end
