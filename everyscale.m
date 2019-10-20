function [MTheta, MPhi] = everyscale(Nx, Ny, MinT, MinP, Mtheta, Mphi,...
    CWx, CWy)
% Функція ппершого масштабування пікселів

% Mtheta, Mphi - координати точок проекції пікселів
% MinT, MinP - масштабні коефіцієнти
% Nx, Ny - кількість точок які можливо відобразити
% CWx, CWy - центра мас для "сусідів" пікселів та межі

% Визначаємо менший із масштабів
if ~isnan(MinT) && ~isnan(MinP)
    MinSc = min(MinT, MinP);
else
    MinSc = 1.0;
end

% Масштабуємо
if MinSc ~= 1.0
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:5
                MTheta{k}(i, j) = (MinSc * (Mtheta{k}(i, j) - CWx{1}(i, j)))...
                    + CWx{1}(i, j);
                MPhi{k}(i, j) = (MinSc * (Mphi{k}(i, j) - CWy{1}(i, j)))...
                    + CWy{1}(i, j);
            end
        end
    end
else
    MTheta = Mtheta;
    MPhi = Mphi;
end
    
end
