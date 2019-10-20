script % Projection_of_pixels_2
% Запуск всього F5
% Запуск виділеного F9

%%%%%%%%
% Дано %
%%%%%%%%
clear; clc;

% Кількість активних елементів розкладу
% в напрямку польоту (по вертикалі, Ox) TDI
pd = 32 + 1;
% поперек польоту (по горизонталі, Oy)
qd = 4096 + 1;

% Номер пікселя

% Розмір чутливого елемента
% довжина [мкм]
vd = 8 * 10 ^ -6;
% ширина [мкм]
wd = 8 * 10 ^ -6;

% Період чутливих елементів
% вздовж [мкм]
Vd = 8 * 10 ^ -6;
% впоперек [мкм]
Wd = 8 * 10 ^ -6;

% Фокусна відстань об'єктива [мм]
f = 2260 * 10 ^ -3;

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [град]
dAp = 0 * (pi / 180.0);
% впоперек [град]
dAq = 0 * (pi / 180.0);

% Висота орбіти  [км]
h = 668 * 10 ^ 3;

% Кути візування
% тангажа [град]
theta = 0 * (pi / 180.0);
% крена [град]
phi = 0 * (pi / 180.0);
% рискання [град]
psi = 0 * (pi / 180.0);

%%%%%%%%%%%%%%%%%%%%%%%%
% Уточнення параметрів %
%%%%%%%%%%%%%%%%%%%%%%%%

% Кількість опорних пікселів для відображення
% вздовж напрямку польоту
Nx0 = 7;

% поперек напрямку польоту
Ny0 = 11;

% Відображення точок
% якщо id = 0 - відображати всі
% якщо id = 1 - відображати тільки вибрані
idx = 1;
idy = 1;
% Проведення масштабування
id = 1;

% Проміжок між пікселями в % (відсотках)
proc = 10;

%%%%%%%%%%%%%%
% РОЗРАХУНОК %
%%%%%%%%%%%%%%

% Корекція кількості точок які необхідно відобразити
Nx = numbertoch(Nx0, pd, idx);
Ny = numbertoch(Ny0, qd, idy);

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [мкм]
dLp = f * tan(dAp);
% впоперек [мкм]
dLq = f * tan(dAq);

clear dAp dAq

% Початкове лінійне значення
Lp0 = dLp - 0.5 * Vd * (pd - 1);
Lq0 = dLq - 0.5 * Wd * (qd - 1);

clear dLp dLq

% Визначаємо нумерацію відображуваних пікселів
Mx = masivtoch(Nx, pd, idx);
My = masivtoch(Ny, qd, idy);

% Корекція кількості точок які можливо відобразити
Nx = length(Mx);
Ny = length(My);

clear pd qd

% Розрахунок лінійних координит центрів пікселів
for i = 1:Nx
    Lpi(i) = Lp0 + (Mx(i) - 1) * Vd;
end

for i = 1:Ny
    Lqi(i) = Lq0 + (My(i) - 1) * Wd;
end

clear i Lp0 Lq0

% Розрахунок кутових координит центрів пікселів
% Всі значення для точки "0" (0, 0)
Wxi{1} = atan(Lpi ./ f);
Wyi{1} = atan(Lqi ./ f);

% Добавляємо центри точок для пікселів (повністю)
WXi{1} = Wxi{1};
WYi{1} = Wyi{1};

clear Lpi Lqi

% Оскільки розмір активної частини пікселя використовується лише для
% енергетичної складової, а просторова роздільна здатність визначається
% всім розміром пікселя, то необхідно розраховувати 2 масиви; також важливо
% замітити, що 2-й масив з розрахунком всього пікселя використовуватиметься
% для алгоритму відображення пікселів (алгоритм ущільнення)

% Кутові величини активної частини (АЧ) пікселів
wx = anglesize(Wxi{1}, vd, f);
wy = anglesize(Wyi{1}, wd, f);

% Кутові величини пікселів (повністю)
Wx = anglesize(Wxi{1}, Vd, f);
Wy = anglesize(Wyi{1}, Wd, f);

% Різниці кутових величин АЧ пікселів відносно їх центрів
dwx = diffangle(Wxi{1}, wx, vd, f);
dwy = diffangle(Wyi{1}, wy, wd, f);

% Різниці кутових величин пікселів (повністю) відносно їх центрів
dWx = diffangle(Wxi{1}, Wx, Vd, f);
dWy = diffangle(Wyi{1}, Wy, Wd, f);

clear vd wd Vd Wd f

% Всі значення для точки 1 (1, -1)
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{2} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 2 (1, 1)
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{3} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 3 (-1, 1)
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{4} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 1 (-1, -1)
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{5} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{5} = koordtoch(Wyi{1}, Wy, my, dWy);

clear mx my wx wy dwx dwy Wx Wy dWx dWy

% В результаті отримали масив даних для точок пікселів
% Врахуємо вплив кута рискання
for k = 1:5
    for i = 1:Nx
       for j = 1:Ny
           % Для АЧ пікселів
           theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
               - tan(Wyi{k}(1, j)) * sin(psi));
           phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
               + tan(Wxi{k}(1, i)) * sin(psi));
           
           % Для всієї частини (ВЧ) пікселів
           Theta0{k}(i, j) = atan(tan(WXi{k}(1, i)) * cos(psi)...
               - tan(WYi{k}(1, j)) * sin(psi));
           Phi0{k}(i, j) = atan(tan(WYi{k}(1, j)) * cos(psi)...
               + tan(WXi{k}(1, i)) * sin(psi));
       end
    end
end

clear psi k i j Wxi Wyi WXi WYi

% Розраховуєм зміну кутових координат в залежності від кутів нахилу
for k = 1:5
    % спочатку по тангажу, а тоді по крену
    [theta1{k}, phi1{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 0);
    [Theta1{k}, Phi1{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 0);
    
    % спочатку по крену, а тоді по тангажу
    [theta2{k}, phi2{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 1);
    [Theta2{k}, Phi2{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 1);
end

clear k theta phi theta0 phi0 Theta0 Phi0

% Розрахунок кутів нахилу колонок і стовбців
% для рядків (поперек польоту)
for i = 1:Nx
    if (tan(theta1{1}(i, 1)) - tan(theta1{1}(i, Ny))) ~= 0
        PhiR1(i) = acot((tan(theta1{1}(i, 1)) - tan(theta1{1}(i, Ny))) \...
            (tan(phi1{1}(i, 1)) - tan(phi1{1}(i, Ny))));
        % В градусах
        PhiR1(i) = 90 - PhiR1(i) * (180 / pi);
    else
        % В градусах
        PhiR1(i) = 90;
    end
    
    if (tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) ~= 0
        PhiR2(i) = acot((tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) \...
            (tan(phi2{1}(i, 1)) - tan(phi2{1}(i, Ny))));
        % В градусах
        PhiR2(i) = 90 - PhiR2(i) * (180 / pi);
    else
        % В градусах
        PhiR2(i) = 90;
    end
end

clear i

% для стовбців (в напрямку польоту)
for j = 1:Ny
    if (tan(theta1{1}(1, j)) - tan(theta1{1}(Nx, j))) ~= 0
        PhiC1(j) = atan((tan(theta1{1}(1, j)) - tan(theta1{1}(Nx, j))) \...
            (tan(phi1{1}(1, j)) - tan(phi1{1}(Nx, j))));
        % В градусах
        PhiC1(j) = PhiC1(j) * (180 / pi);
    else
        % В градусах
        PhiC1(j) = 90;
    end
    
    if (tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) ~= 0
        PhiC2(j) = atan((tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) \...
            (tan(phi2{1}(1, j)) - tan(phi2{1}(Nx, j))));
        % В градусах
        PhiC2(j) = PhiC2(j) * (180 / pi);
    else
        % В градусах
        PhiC2(j) = 90;
    end
end

clear j

% Розраховуємо лінійні координати точок
for k = 1:5
    % Для АЧ пікселів
    Mtheta1{k} = h .* tan(theta1{k});
    Mtheta2{k} = h .* tan(theta2{k});
    
    Mphi1{k} = h .* tan(phi1{k});
    Mphi2{k} = h .* tan(phi2{k});
    
    % Для ВЧ пікселів
    MTheta1{k} = h .* tan(Theta1{k});
    MTheta2{k} = h .* tan(Theta2{k});
    
    MPhi1{k} = h .* tan(Phi1{k});
    MPhi2{k} = h .* tan(Phi2{k});
end

clear h k Theta1 Phi1 Theta2 Phi2 theta1 phi1 theta2 phi2 


% Відображення тільки необхідних пікселів, а не всього масиву
if (id == 1) && (idx == 1) && (idy == 1) && (Nx0 > 1) && (Ny0 > 1)
    % АЛГОРИТМ УЩІЛЬНЕННЯ
    % Розрахунок нумерацій точок з врахуванням сусідів
    Mx = sysidtoch(Mx);
    My = sysidtoch(My);

    % Знаходимо центр мас для пікселів (center weight - CW)
    [CWx1, CWy1] = centerweight(Mx, My, MTheta1, MPhi1, Nx, Ny);
    [CWx2, CWy2] = centerweight(Mx, My, MTheta2, MPhi2, Nx, Ny);

    % Утворюємо матрицю яка дозволить визначити мінамільні коефіцієнти
    [MinT1, MinP1] = minvidmizpix(Mx, My, CWx1, CWy1, Nx, Ny, proc);
    [MinT2, MinP2] = minvidmizpix(Mx, My, CWx2, CWy2, Nx, Ny, proc);
    
    % Масштабуємо по найменшому коефіцієнтові кожен піксель
    [Mtheta1, Mphi1] = everyscale(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    % Той коефіцієнт який не дорівнює 1 масштабуватиме, але виключно або
    % рядки або стовбці, притому необхідно визначити нові центри мас
    [CWx1, CWy1] = centerweight2(MTheta1, MPhi1, MinT1, MinP1,...
        Nx, Ny, CWx1, CWy1);
    [CWx2, CWy2] = centerweight2(MTheta2, MPhi2, MinT2, MinP2,...
        Nx, Ny, CWx2, CWy2);
    
    % Необхідно переаналізувати отримані дані, для повторного масштабування
    % по іншій із осей
    % Скорегуємо коефіцієнти після масштабування
    Min1 = min(MinT1, MinP1);
    MinT1 = MinT1 / Min1;
    MinP1 = MinP1 / Min1;
    Min2 = min(MinT2, MinP2);
    MinT2 = MinT2 / Min2;
    MinP2 = MinP2 / Min2;
    
    clear Min1 Min2
    
    % Масштабуємо по найменшому коефіцієнтові кожен піксель
    [Mtheta1, Mphi1] = everyscale2(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale2(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    clear MinT1 MinP1 MinT2 MinP2 %CWx1 CWy1 CWx2 CWy2
end

clear idx idy Nx0 Ny0 Mx My

% Малюємо зображення проекцій
figure(1);
clf;
hold on;
grid off;   % включити сітку "on", виключити "off"
xlabel('Поперек польоту [м]');
ylabel('Вздовж польоту [м]');
title('Проекція матриці пікселів (тангаж - крен)');
% plot(0,0, 'k^');

% Визначаємо максимум і мінімум для однакового масштабу по осям
xmin = min(Mtheta1{1}(:));
xmax = max(Mtheta1{1}(:));
ymin = min(Mphi1{1}(:));
ymax = max(Mphi1{1}(:));

for k = 2:5
    ymin = min(min(Mphi1{k}(:)), ymin);
    ymax = max(max(Mphi1{k}(:)), ymax);
    xmin = min(min(Mtheta1{k}(:)), xmin);
    xmax = max(max(Mtheta1{k}(:)), xmax);
end

for k = 1:5
    ymin = min(min(Mphi2{k}(:)), ymin);
    ymax = max(max(Mphi2{k}(:)), ymax);
    xmin = min(min(Mtheta2{k}(:)), xmin);
    xmax = max(max(Mtheta2{k}(:)), xmax);
end

% Масштабуємо координатні осі
ylen = ymax - ymin;
xlen = xmax - xmin;

ycenter = ymin + 0.5 * ylen;
xcenter = xmin + 0.5 * xlen;

clear ymax ymin xmax xmin

len = max(xlen, ylen) * (1 + 0.5 * proc / 100);

clear xlen ylen proc

clear xlen ylen

ylim(xcenter - 0.5 .* [len -len]);
xlim(ycenter - 0.5 .* [len -len]);

%  Відображаємо проекції пікселів
for i = 1:Nx
   for j = 1:Ny
       plot(Mphi1{1}(i, j), Mtheta1{1}(i, j), 'r+');
       plot([Mphi1{2}(i, j), Mphi1{3}(i, j),...
           Mphi1{4}(i, j), Mphi1{5}(i, j), Mphi1{2}(i, j)],...
           [Mtheta1{2}(i, j), Mtheta1{3}(i, j), Mtheta1{4}(i, j),...
           Mtheta1{5}(i, j), Mtheta1{2}(i, j)], '-k');
%        if id == 1
%            plot(CWy1{1}(i, j), CWx1{1}(i, j), 'b+');
%        end
   end
end

% Перевірка чи лежать центри пікселів на одній лінії
% plot([Mphi1{1}(Nx, 1), Mphi1{1}(Nx, Ny),],...
%     [Mtheta1{1}(Nx, 1), Mtheta1{1}(Nx, Ny)], '-g');
% plot([Mphi1{1}(1, 1), Mphi1{1}(Nx, 1),],...
%     [Mtheta1{1}(1, 1), Mtheta1{1}(Nx, 1)], '-g');
% 
% plot([Mphi1{1}(1, 1), Mphi1{1}(1, Ny),],...
%     [Mtheta1{1}(1, 1), Mtheta1{1}(1, Ny)], '-g');

figure(2);
clf;
hold on;
grid off;   % включити сітку "on", виключити "off"
xlabel('Поперек польоту [м]');
ylabel('Вздовж польоту [м]');
title('Проекція матриці пікселів (крен - тангаж)');

% ylim([ymin ymax]);
% xlim([xmin xmax]);

ylim(xcenter - 0.5 .* [len -len]);
xlim(ycenter - 0.5 .* [len -len]);

clear xcenter ycenter len

%  Відображаємо проекції пікселів
for i = 1:Nx
   for j = 1:Ny
       plot(Mphi2{1}(i, j), Mtheta2{1}(i, j), 'r+');
       plot([Mphi2{2}(i, j), Mphi2{3}(i, j),...
           Mphi2{4}(i, j), Mphi2{5}(i, j), Mphi2{2}(i, j)],...
           [Mtheta2{2}(i, j), Mtheta2{3}(i, j), Mtheta2{4}(i, j),...
           Mtheta2{5}(i, j), Mtheta2{2}(i, j)], '-k');
%        if id == 1
%            plot(CWy2{1}(i, j), CWx2{1}(i, j), 'b+');
%        end
   end
end

% Перевірка чи лежать центри пікселів на одній лінії
% plot([Mphi2{1}(Nx, 1), Mphi2{1}(Nx, Ny),],...
%     [Mtheta2{1}(Nx, 1), Mtheta2{1}(Nx, Ny)], '-g');
% plot([Mphi2{1}(1, 1), Mphi2{1}(Nx, 1),],...
%     [Mtheta2{1}(1, 1), Mtheta2{1}(Nx, 1)], '-g');
% 
% plot([Mphi2{1}(1, 1), Mphi2{1}(1, Ny),],...
%     [Mtheta2{1}(1, 1), Mtheta2{1}(1, Ny)], '-g');

clear k i j pd qd id CWx1 CWy1 CWx2 CWy2 Mphi1 Mphi2 Mtheta1 Mtheta2

% Розразунок точок які визначатимуть роздільну здатність системи
% нахил спочатку по тангажу, а тоді по крену
[Ltheta1, Lphi1] = razreshtoch(MTheta1, MPhi1, Nx, Ny);
[Ltheta2, Lphi2] = razreshtoch(MTheta2, MPhi2, Nx, Ny);

clear MTheta1 MPhi1 MTheta2 MPhi2 Nx Ny

% Розрахунок частоти [мм ^ -1]
NuTheta1 = 1000 ./ Ltheta1;
NuLphi1 = 1000 ./ Lphi1;
NuTheta2 = 1000 ./ Ltheta2;
NuLphi2 = 1000 ./ Lphi2;

% figure(1);
% hold on;
% %  Відображаємо проекції пікселів
% for i = 1:Nx
%    for j = 1:Ny
%        plot(Mphi1{1}(i, j), Mtheta1{1}(i, j), 'bo');
%        plot(Mphi1{2}(i, j), Mtheta1{2}(i, j), 'bo');
%        plot(Mphi1{3}(i, j), Mtheta1{3}(i, j), 'bo');
%        plot(Mphi1{4}(i, j), Mtheta1{4}(i, j), 'bo');
%    end
% end
% 
% clear i j
% 
% figure(2);
% hold on;
% %  Відображаємо проекції пікселів
% for i = 1:Nx
%    for j = 1:Ny
%        plot(Mphi2{1}(i, j), Mtheta2{1}(i, j), 'bo');
%        plot(Mphi2{2}(i, j), Mtheta2{2}(i, j), 'bo');
%        plot(Mphi2{3}(i, j), Mtheta2{3}(i, j), 'bo');
%        plot(Mphi2{4}(i, j), Mtheta2{4}(i, j), 'bo');
%    end
% end
% 
% clear i j Nx Ny Mtheta1 Mphi1 Mtheta2 Mphi2

clc;
