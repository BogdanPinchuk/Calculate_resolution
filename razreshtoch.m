function [Ltheta, Lphi] = razreshtoch(MTheta, MPhi, Nx, Ny)
%[Ltheta, Lphi, Mtheta, Mphi] = razreshtoch(MTheta, MPhi, Nx, Ny)
% Функція розразуноку точок які визначатимуть роздільну здатність системи
% для нахил спочатку по тангажу, а тоді по крену і
% для нахил спочатку по крену, а тоді по тангажу
% Та визначення роздільної здатності системи

% Mtheta, Mphi - шукані точки
% Ltheta, Lphi - роздільна здатність системи (грубо вздовж і в поперек польоту)
% MTheta, MPhi - дані масиви точок,
% між якими треба найти відстань

% Знаходимо перетин точок
for i = 1:Nx
    for j = 1:Ny
        [Mtheta{1}(i, j), Mphi{1}(i, j)] = peretintoch(MTheta{2}(i, j),...
            MPhi{2}(i, j), MTheta{3}(i, j), MPhi{3}(i, j), MTheta{1}(1, j),...
            MPhi{1}(1, j), MTheta{1}(Nx, j), MPhi{1}(Nx, j));
        
        [Mtheta{3}(i, j), Mphi{3}(i, j)] = peretintoch(MTheta{4}(i, j),...
            MPhi{4}(i, j), MTheta{5}(i, j), MPhi{5}(i, j), MTheta{1}(1, j),...
            MPhi{1}(1, j), MTheta{1}(Nx, j), MPhi{1}(Nx, j));
        
        [Mtheta{2}(i, j), Mphi{2}(i, j)] = peretintoch(MTheta{3}(i, j),...
            MPhi{3}(i, j), MTheta{4}(i, j), MPhi{4}(i, j), MTheta{1}(i, 1),...
            MPhi{1}(i, 1), MTheta{1}(i, Ny), MPhi{1}(i, Ny));
        
        [Mtheta{4}(i, j), Mphi{4}(i, j)] = peretintoch(MTheta{2}(i, j),...
            MPhi{2}(i, j), MTheta{5}(i, j), MPhi{5}(i, j), MTheta{1}(i, 1),...
            MPhi{1}(i, 1), MTheta{1}(i, Ny), MPhi{1}(i, Ny));
    end
end

% Знаходмо відстані між точками
Ltheta = vipstanmizpix(Mtheta{1}, Mtheta{3}, Mphi{1}, Mphi{3});
Lphi = vipstanmizpix(Mtheta{2}, Mtheta{4}, Mphi{2}, Mphi{4});

end
