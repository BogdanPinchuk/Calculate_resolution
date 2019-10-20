function [MTheta, MPhi] = everyscale3(MaxSc, Mtheta, Mphi)
% Функція кінечного масштабування пікселів

% MaxSc - Результуючий максимальний коеф.
% Mtheta - координати точок проекції пікселів (в напрямку польоту)
% Mphi - координати точок проекції пікселів (поперек польоту)

% Визначаємо розмірність
Nx = size(Mtheta{1}, 1);
Ny = size(Mtheta{1}, 2);

% Найдемо границі для центрування зображення
for k = 1:5
    minT1(k) = min(min(Mtheta{k}));
    maxT1(k) = max(max(Mtheta{k}));
    minP1(k) = min(min(Mphi{k}));
    maxP1(k) = max(max(Mphi{k}));
end

minT = min(minT1);
maxT = max(maxT1);
minP = min(minP1);
maxP = max(maxP1);

% Знаходимо центри
centerT = minT + 0.5 * (maxT - minT);
centerP = minP + 0.5 * (maxP - minP);

% Масштабуємо
for i = 1:Nx
    for j = 1:Ny
        for k = 1:5
            % Зміщуємо - центруємо точки
            MTheta{k}(i, j) = Mtheta{k}(i, j) - centerT;
            MPhi{k}(i, j) = Mphi{k}(i, j) - centerP;
            
            % Масштабуємо
            MTheta{k}(i, j) = MTheta{k}(i, j) / MaxSc;
            MPhi{k}(i, j) = MPhi{k}(i, j) / MaxSc;
        end
    end
end


end
