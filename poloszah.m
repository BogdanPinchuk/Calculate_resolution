function [LzahT, LzahP, CzahT, CzahP] = poloszah(Mtheta, Mphi, Nx, Ny)
% Функція розрахунку номерів точок для визначення полоси захвату

% Mtheta, Mphi - лінійні координати точок
% LzahT, LzahP - полоси захвату вздовж і поперек польту

Lmin = min([Mphi{2}(Nx, 1), Mphi{3}(Nx, Ny), Mphi{4}(1, Ny), Mphi{5}(1, 1)]);
Lmax = max([Mphi{2}(Nx, 1), Mphi{3}(Nx, Ny), Mphi{4}(1, Ny), Mphi{5}(1, 1)]);

LzahP = Lmax - Lmin;
CzahP = Lmin + 0.5 .* LzahP;

Lmin = min([Mtheta{2}(Nx, 1), Mtheta{3}(Nx, Ny), Mtheta{4}(1, Ny), Mtheta{5}(1, 1)]);
Lmax = max([Mtheta{2}(Nx, 1), Mtheta{3}(Nx, Ny), Mtheta{4}(1, Ny), Mtheta{5}(1, 1)]);

LzahT = Lmax - Lmin;
CzahT = Lmin + 0.5 .* LzahT;

end
