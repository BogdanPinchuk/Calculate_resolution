script % Projection_of_pixels_7
% Запуск всього F5
% Запуск виділеного F9

%%%%%%%%
% Дано %
%%%%%%%%
clear; clc;

% Кількість активних елементів розкладу
% в напрямку польоту (по вертикалі, Ox) TDI
pd = 128;

% поперек польоту (по горизонталі, Oy)
qd = 24001;

% Номер пікселя

% Розмір чутливого елемента
% довжина [мкм]
vd = 7 * 10 ^ -6;  %8

% ширина [мкм]
wd = 7 * 10 ^ -6;  %8

% Період чутливих елементів
% вздовж [мкм]
Vd = 7 * 10 ^ -6;  %8
% впоперек [мкм]
Wd = 7 * 10 ^ -6;  %8

% Фокусна відстань об'єктива [мм]
f = 6680  * 10 ^ -3; %2260

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [град]
dAp = 0 * (pi / 180.0);
% впоперек [град]
dAq = 0 * (pi / 180.0);

% Висота орбіти [км]
h = 668 * 10 ^ 3;

% Кути візування
% тангажа [град]
Theta = 35 * (pi / 180.0);
% крена [град]
Phi = 35 * (pi / 180.0);
% рискання [град]
psi = 0 * (pi / 180.0);

% Параметри для аналізу відхилення точності формул
% Кількість точок
Ntk = 8;

% Широта (-90°...+90°) на якій знаходиться КА [град]
gama = 45 * (pi / 180.0);

% Радіуси Землі, тобто відстань від точки поверхні Землі до її центру [км]
Rz = 6371.032 * 10 ^ 3;
Rzmin = 6356.777 * 10 ^ 3;
Rzmax = 6378.160 * 10 ^ 3;

%%%%%%%%%%%%%%%%%%%%%%%%
% Уточнення параметрів %
%%%%%%%%%%%%%%%%%%%%%%%%

% Кількість опорних пікселів для відображення
% вздовж напрямку польоту
Nx0 = 3;

% поперек напрямку польоту
Ny0 = 3;

% Відображення точок
% якщо id = 0 - відображати всі
% якщо id = 1 - відображати тільки вибрані
idx = 1;
idy = 1;
% Проведення масштабування (0 - ні, 1 - так)
id = 1;

% Проміжок між пікселями в % (відсотках)
proc = 10;

% Відображення графіків
graph = false;

% Варіанти розрахунку 1 (false) або 
% більше пікселів (true)
Npix = true;

% Розрахунок або все true колонки чи рядка false
CR = true;

% Розрахунок колонки true, або рядка false
HV = false;

% Номер колонки або рядка
numHV = 20000;
 
 % Нумерація вибраного пікселя (зліва на право, знизу до верху)
pixNx = 20;
pixNy = 20000;
 
%%%%%%%%%%%%%%
% РОЗРАХУНОК %
%%%%%%%%%%%%%%

% Корекція кількості точок які необхідно відобразити
if Npix % розрахунок 1 пікселя або всіх
    if CR % розрахунок всього або окремої колонки
        Nx = numbertoch(Nx0, pd, idx);
        Ny = numbertoch(Ny0, qd, idy);
    else
        if HV % розрахунок колонки
            Nx = numbertoch(Nx0, pd, idx);
            Ny = 1;
        else % розрахунок рядка
            Nx = 1;
            Ny = numbertoch(Ny0, qd, idy);
        end
    end
else
    Nx = 1;
    Ny = 1;
end

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [мкм]
dLp = f * tan(dAp);
% впоперек [мкм]
dLq = f * tan(dAq);

% clear dAp dAq

% Початкове лінійне значення
Lp0 = 0.5 * Vd * (1 - pd) - dLp;
Lq0 = 0.5 * Wd * (1 - qd) - dLq;

clear dLp dLq

% Визначаємо нумерацію відображуваних пікселів
if Npix % розрахунок 1 пікселя або всіх
    if CR % розрахунок всього або окремої колонки
        Mx = masivtoch(Nx, pd, idx);
        My = masivtoch(Ny, qd, idy);
    else
        if HV % розрахунок колонки
            Mx = masivtoch(Nx, pd, idx);
            My = numHV;
        else % розрахунок рядка
            Mx = numHV;
            My = masivtoch(Ny, qd, idy);
        end
    end
else
    Mx = pixNx;
    My = pixNy;
end

% Корекція кількості точок які можливо відобразити
Nx = length(Mx);
Ny = length(My);

clear CR Npix pixNx pixNy HV numHV

% Розрахунок лінійних координит центрів пікселів
for i = 1:Nx
%     Lpi(i) = Lp0 + (Mx(i) - 1) * Vd;
    
    % Перевернемо масив, щоб скорегувати дані відносно розташування на ЗП
    Lpi(i) = Lp0 + (Mx(Nx + 1 - i) - 1) * Vd;
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
% зауважити, що 2-й масив з розрахунком всього пікселя використовуватиметься
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

% clear vd wd

% Розрахунок активної частини АЧ (Wxi, Wyi)
% Розрахунок всієї частини ВС (WXi, Wyi)
% Розрахунок для роздільної здатності (wxi, wyi)
% Всі значення для точки 1 (1, -1) та (1, 0) 5
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
wxi{1} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{2} = koordtoch(Wyi{1}, Wy, my, dWy);
my = 0;
wyi{1} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 2 (1, 1) та (0, 1) 6
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
mx = 0;
wxi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{3} = koordtoch(Wyi{1}, Wy, my, dWy);
wyi{2} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 3 (-1, 1) та (-1, 0) 7
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
wxi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{4} = koordtoch(Wyi{1}, Wy, my, dWy);
my = 0;
wyi{3} = koordtoch(Wyi{1}, Wy, my, dWy);

% Всі значення для точки 1 (-1, -1) та (0, -1) 8
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{5} = koordtoch(Wxi{1}, Wx, mx, dWx);
mx = 0;
wxi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{5} = koordtoch(Wyi{1}, Wy, my, dWy);
wyi{4} = koordtoch(Wyi{1}, Wy, my, dWy);

clear mx my wx wy dwx dwy dWx dWy Wx Wy

% В результаті отримали масив даних для точок пікселів
% Врахуємо вплив кута рискання
for k = 1:5
    for i = 1:Nx
       for j = 1:Ny
           % Для АЧ пікселів
%            theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
%                - tan(Wyi{k}(1, j)) * sin(psi));
%            phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
%                + tan(Wxi{k}(1, i)) * sin(psi));
%            
%            % Для всієї частини (ВЧ) пікселів
%            Theta0{k}(i, j) = atan(tan(WXi{k}(1, i)) * cos(psi)...
%                - tan(WYi{k}(1, j)) * sin(psi));
%            Phi0{k}(i, j) = atan(tan(WYi{k}(1, j)) * cos(psi)...
%                + tan(WXi{k}(1, i)) * sin(psi));
           
           % Для роздільної здатності
           if k < 5
               Rtheta0{k}(i, j) = atan(tan(wxi{k}(1, i)) * cos(psi)...
                   - tan(wyi{k}(1, j)) * sin(psi));
               Rphi0{k}(i, j) = atan(tan(wyi{k}(1, j)) * cos(psi)...
                   + tan(wxi{k}(1, i)) * sin(psi));
           end
       end
    end
end

clear k i j Wxi Wyi wxi wyi %WXi WYi psi

% Перед тим як розраховувати вплив кутів нахил, розрахуємо вплив кривизни Землі 
% на кути нахилу (тангаж і крена)
            
%////////////////////
%// Кривизна Землі //
%////////////////////

% Розрахунок роздільної здатності враховуючи кривизну Землі
% Розраховуємо відстань від центра до точкин на поверхні Землі
Rt = radzem(Rzmin, Rzmax, gama);

% Розраховуємо радіус кривизни Землі
Rk = radkrivzem(Rzmin, Rzmax, gama);

% clear gama

% Розраховуємо похибку на висоту
dh = Rt - Rz;

clear Rt Rzmin Rzmax Rz

% Корегуємо висоту
H = h + dh;

clear dh %h

%%%%%%%%%%%%%%%%%%%%%%%%%
% Початок
% Циклу розрахунку масиву значень відхилень за тангажем і креном
%%%%%%%%%%%%%%%%%%%%%%%%%

% крок тангажа
thetaE = Theta / (Ntk - 1);
% крок крена
phiE = Phi / (Ntk - 1);

for it = 1:Ntk
    for ik = 1:Ntk
        theta = (it - 1) * thetaE;
        phi = (ik - 1) * phiE;
        
        % Корегуємо висоту із врахуванням впливу на кути відхилень за рахунок кривизни Землі
        % тобто це висота до площини сканування (площина дотичної до точки відхилення КА)
        Hrz = heighttoploschin(H, Rk, theta, phi);

        % Кути між двома площинами, в перерізі по площині в напрямку (0°)
        % та поперек польту (90°) і основною площиною
        deltaTheta = kytmizploschinami(0.0, theta, phi, H, Rk);
        deltaPhi = kytmizploschinami(0.5 * pi, theta, phi, H, Rk);

        % Корегуємо кути візування згідно із врахуванням кривизни Землі
        theta_0 = theta + deltaTheta;
        phi_0 = phi + deltaPhi;

        clear deltaPhi deltaTheta

        % Розраховуємо еквівалентні координати центра проекції МПВ
        [l_theta, l_phi] = centerMPV(H, Rk, theta, phi, theta_0, phi_0);

%         clear H Rk

        % Різниця довжин між розрахованим (із зміненою висотою) і реальною довжиною по дузі
        deltaL_t = Hrz * tan(theta_0) - l_theta;
        deltaL_p = Hrz * tan(phi_0) - l_phi;

        clear l_phi l_theta

        % Розраховуєм зміну кутових координат в залежності від кутів нахилу
        for k = 1:5
            if k < 5
                % спочатку по тангажу, а тоді по крену
        %         [Rtheta1{k}, Rphi1{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
                % спочатку по крену, а тоді по тангажу
        %         [Rtheta2{k}, Rphi2{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);

                % додаткові формули для локального повороту
                % спочатку по тангажу, а тоді по крену
                [Rtheta3{k}, Rphi3{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
                % спочатку по крену, а тоді по тангажу
%                 [Rtheta4{k}, Rphi4{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);
            end
        end

        clear k %theta0 phi0 Rtheta0 Rphi0
            %phi_0 theta_0 Theta0 Phi0 theta phi

        % Розраховуємо лінійні координати точок
        for k = 1:5
            % Для РЗ роздільної здатності
            if k < 5
        %         MRtheta1{k} = Hrz .* tan(Rtheta1{k}) + deltaL_t;
        %         MRtheta2{k} = Hrz .* tan(Rtheta2{k}) + deltaL_t;
        %         
        %         MRphi1{k} = Hrz .* tan(Rphi1{k}) + deltaL_p;
        %         MRphi2{k} = Hrz .* tan(Rphi2{k}) + deltaL_p;

                % додаткові формули для локального повороту
                MRtheta3{k} = Hrz .* tan(Rtheta3{k}) + deltaL_t;
%                 MRtheta4{k} = Hrz .* tan(Rtheta4{k}) + deltaL_t;

                MRphi3{k} = Hrz .* tan(Rphi3{k}) + deltaL_p;
%                 MRphi4{k} = Hrz .* tan(Rphi4{k}) + deltaL_p;
            end
        end

        clear k Theta1 Phi1 Theta2 Phi2 theta1 phi1 theta2 phi2...
            Theta3 Phi3 Theta4 Phi4 theta3 phi3 theta4 phi4...
            deltaL_p deltaL_t ...
            Rphi1 Rphi2 Rphi3 Rphi4 Rtheta1 Rtheta2 Rtheta3 Rtheta4

        % Розразунок точок які визначатимуть роздільну здатність системи (на площині)
        % спочатку по тангажу, а тоді по крену
        % Ltheta1 = vipstanmizpix(MRtheta1{1}, MRtheta1{3}, MRphi1{1}, MRphi1{3});
        % Lphi1 = vipstanmizpix(MRtheta1{2}, MRtheta1{4}, MRphi1{2}, MRphi1{4});

        Ltheta3 = vipstanmizpix(MRtheta3{1}, MRtheta3{3}, MRphi3{1}, MRphi3{3});
        Lphi3 = vipstanmizpix(MRtheta3{2}, MRtheta3{4}, MRphi3{2}, MRphi3{4});

        % спочатку по крену, а тоді по тангажу
        % Ltheta2 = vipstanmizpix(MRtheta2{1}, MRtheta2{3}, MRphi2{1}, MRphi2{3});
        % Lphi2 = vipstanmizpix(MRtheta2{2}, MRtheta2{4}, MRphi2{2}, MRphi2{4});

%         Ltheta4 = vipstanmizpix(MRtheta4{1}, MRtheta4{3}, MRphi4{1}, MRphi4{3});
%         Lphi4 = vipstanmizpix(MRtheta4{2}, MRtheta4{4}, MRphi4{2}, MRphi4{4});
        
        clear MTheta3 MPhi3 MTheta4 MPhi4 MRtheta3 MRphi3 MRtheta4 MRphi4...
            MRphi1 MRphi2 MRtheta1 MRtheta2

        % Розрахунок просторового розділення за спрощеними фонрмулами Тягур В.М.
        % в напрямку польоту
        A = Vd .* Hrz ./ f;
        % поперек польоту
        P = Wd .* Hrz ./ f;

        clear Hrz %Vd Wd f

        % Додаткові кути тетраедра при відхиленні
        theta_2 = atan(tan(phi_0) .* cos(theta_0));
        phi_2 = atan(tan(theta_0) .* cos(phi_0));

        for i = 1:Nx
           for j = 1:Ny
               % Для роздільної здатності

               %%%%%%%%%%%%%%%%%%%%%
               % Система тангаж-крен (лок. пов.)
               % Система крен-тангаж (глоб. пов.)
               %%%%%%%%%%%%%%%%%%%%%

               % в напрямку польоту
        %        A3(i, j) = A ./ (cos(theta_0 + Theta0{1}(i, j)) .^ 2 .* ...
        %            cos(theta_2 + Phi0{1}(i, j)));
               A3(i, j) = A ./ (cos(theta_0 + WXi{1}(i)) .^ 2 .* ...
                   cos(theta_2 + WYi{1}(j)));
               % поперек польоту
        %        P3(i, j) = P ./ (cos(theta_0 + Theta0{1}(i, j)) .* ...
        %            cos(theta_2 + Phi0{1}(i, j)) .^ 2);
               P3(i, j) = P ./ (cos(theta_0 + WXi{1}(i)) .* ...
                   cos(theta_2 + WYi{1}(j)) .^ 2);

               %%%%%%%%%%%%%%%%%%%%%
               % Система крен-тангаж (лок. пов.)
               % Система тангаж-крен (глоб. пов.)
               %%%%%%%%%%%%%%%%%%%%%

               % в напрямку польоту
        %        A4(i, j) = A ./ (cos(phi_2 + Theta0{1}(i, j)) .^2 .* ...
        %            cos(phi_0 + Phi0{1}(i, j)));
%                A4(i, j) = A ./ (cos(phi_2 + WXi{1}(i)) .^2 .* ...
%                    cos(phi_0 + WYi{1}(j)));
               % поперек польоту
        %        P4(i, j) = P ./ (cos(phi_2 + Theta0{1}(i, j)) .* ...
        %            cos(phi_0 + Phi0{1}(i, j)) .^2);
%                P4(i, j) = P ./ (cos(phi_2 + WXi{1}(i)) .* ...
%                    cos(phi_0 + WYi{1}(j)) .^2);
           end
        end

        clear theta_2 phi_2 theta_0 phi_0 A P Theta0 Phi0 %WXi WYi
        clear idx idy Nx0 Ny0 Mx My
        
        % Абсолютне відхилення        
        dAbs3t = abs(Ltheta3 - A3);
        dAbs3k = abs(Lphi3 - P3);
%         dAbs4t = abs(Ltheta4 - A4);
%         dAbs4k = abs(Lphi4 - P4);
        
        % Відносне відхилення
        dVid3t = 100 .* abs(1.0 - (A3 ./ Ltheta3));
        dVid3k = 100 .* abs(1.0 - (P3 ./ Lphi3));
%         dVid4t = 100 .* (1.0 - abs(A4 ./ Ltheta4));
%         dVid4k = 100 .* (1.0 - abs(P4 ./ Lphi4));
        
        % Заносимо в глобальну змінну для виведення
        deltaAbs3t{it, ik} = dAbs3t;
        deltaAbs3k{it, ik} = dAbs3k;
%         deltaAbs4t{it, ik} = dAbs4t;
%         deltaAbs4k{it, ik} = dAbs4k;
        
        deltaVid3t{it, ik} = dVid3t;
        deltaVid3k{it, ik} = dVid3k;
%         deltaVid4t{it, ik} = dVid4t;
%         deltaVid4k{it, ik} = dVid4k;
        
        clear A3 A4 P3 P4 Lphi3 Lphi4 Ltheta3 Ltheta4...
         dAbs3k dAbs3t dAbs4k dAbs4t dVid3k dVid3t dVid4k dVid4t
    end
end

clear H Rk WXi WYi Rphi0 Rtheta0 it ik

% Переформатування даних для відображення
for i = 1:Nx
   for j = 1:Ny
       for it = 1:Ntk
           theta = round((it - 1) * thetaE * (180.0 / pi), 5);
           
           Dabs3t{i, j}(it + 1, 1) = theta;
           Dabs3k{i, j}(it + 1, 1) = theta;
           Dvid3t{i, j}(it + 1, 1) = theta;
           Dvid3k{i, j}(it + 1, 1) = theta;
       end
       
       for ik = 1:Ntk
           phi = round((ik - 1) * phiE * (180.0 / pi), 5);
           
           Dabs3t{i, j}(1, ik + 1) = phi;
           Dabs3k{i, j}(1, ik + 1) = phi;
           Dvid3t{i, j}(1, ik + 1) = phi;
           Dvid3k{i, j}(1, ik + 1) = phi;
       end
       
       for it = 1:Ntk
           for ik = 1:Ntk
               Dabs3t{i, j}(it + 1, ik + 1) = deltaAbs3t{it, ik}(i, j);
               Dabs3k{i, j}(it + 1, ik + 1) = deltaAbs3k{it, ik}(i, j);
               Dvid3t{i, j}(it + 1, ik + 1) = deltaVid3t{it, ik}(i, j);
               Dvid3k{i, j}(it + 1, ik + 1) = deltaVid3k{it, ik}(i, j);
           end
       end
   end
end

% clear phiE thetaE
clear it ik %deltaAbs3t deltaAbs3k deltaVid3t deltaVid3k

%%%%%%%%%%%%%%%%%%%%%%%%%
% Кінець
% Циклу розрахунку масиву значень відхилень за тангажем і креном
%%%%%%%%%%%%%%%%%%%%%%%%%

clear k i j id CWx1 CWy1 CWx2 CWy2 Mphi1 Mphi2 Mtheta1 Mtheta2...
    CWx3 CWy3 CWx4 CWy4 Mphi3 Mphi4 Mtheta3 Mtheta4...
    len xcenter ycenter
clc;
clear ans graph proc

% Експеремент
% назва файла
filename = [strrep(['Звіт по розрахунку відхилень ', datestr(now)],...
    ':', '.'), '.xlsx'];

% filename = 'Test.xlsx';

% формуємо для першої сторінки звіт вхідних параметрів
% Змінна для звіту
zvit = {};
zvit{size(zvit, 1) + 1, 2} = strrep(['Звіт по розрахунку відхилень ',...
    datestr(now)], ':', '.');
zvit{size(zvit, 1) + 1, 2} = 'Порядок відхилення для ЛСК - тангаж-крен';
zvit{size(zvit, 1) + 1, 1} = [];

zvit{size(zvit, 1) + 1, 1} = pd;
zvit{size(zvit, 1), 2} = 'Кількість активних елементів розкладу по вертикалі';
zvit{size(zvit, 1) + 1, 1} = qd;
zvit{size(zvit, 1), 2} = 'Кількість активних елементів розкладу по горизонталі';
zvit{size(zvit, 1) + 1, 1} = vd * 10^6;
zvit{size(zvit, 1), 2} = 'Розмір чутливого елемента, довжина [мкм]';
zvit{size(zvit, 1) + 1, 1} = wd * 10^6;
zvit{size(zvit, 1), 2} = 'Розмір чутливого елемента, ширина [мкм]';
zvit{size(zvit, 1) + 1, 1} = Vd * 10^6;
zvit{size(zvit, 1), 2} = 'Період чутливих елементів, вздовж [мкм]';
zvit{size(zvit, 1) + 1, 1} = Wd * 10^6;
zvit{size(zvit, 1), 2} = 'Розмір чутливого елемента, впоперек [мкм]';
zvit{size(zvit, 1) + 1, 1} = f * 10^3;
zvit{size(zvit, 1), 2} = 'Фокусна відстань об`єктива [мм]';
zvit{size(zvit, 1) + 1, 1} = dAp * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Зміщення вздовж осей координат, вздовж [град]';
zvit{size(zvit, 1) + 1, 1} = dAq * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Зміщення вздовж осей координат, впоперек [град]';
zvit{size(zvit, 1) + 1, 1} = h * 10^(-3);
zvit{size(zvit, 1), 2} = 'Висота орбіти [км]';
zvit{size(zvit, 1) + 1, 1} = Theta * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Кут візування, тангаж [град]';
zvit{size(zvit, 1) + 1, 1} = Phi * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Кут візування, крен [град]';
zvit{size(zvit, 1) + 1, 1} = psi * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Кут візування, рискання [град]';
zvit{size(zvit, 1) + 1, 1} = gama * (180.0 / pi);
zvit{size(zvit, 1), 2} = 'Широта Землі на якій знаходиться КА [град]';
zvit{size(zvit, 1) + 1, 1} = [];

zvit{size(zvit, 1) + 1, 2} = 'Порядок пікселів відповідно до проекції на землі';
zvit{size(zvit, 1) + 1, 2} = 'В напрямку польоту 1..Np знизу до верху';
zvit{size(zvit, 1) + 1, 2} = 'Поперек польоту 1..Nq зліва на право';

clear pd qd vd wd Vd Wd f dAp dAq h theta phi psi gama

% Записуємо у файл
xlswrite(filename, zvit);

% Присвоюємо відповідно до точки і зберінаємо кожен пікселя на окремому
% листку
for i = 1:Nx
   for j = 1:Ny
       % Чистимо для подальших вивидень значень
       zvit = {};
       
       y = 1;
       x = 1;
       
       % Відносні значення в напрямку польоту
       zvit{x, y} =...
           'Відносне відхилення в напрямку польоту, [%]';

       x = 2;
       
       % Записуємо дані
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dvid3t{i, j}(it, ik);
          end
       end
       
       y = Ntk + 1 + 2;
       x = 1;
       
       % Відносні значення в поперек польоту
       zvit{x, y} =...
           'Відносне відхилення в поперек польоту, [%]';
       
       x = 2;
       
       % Записуємо дані
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dvid3k{i, j}(it, ik);
          end
       end
       
       y = 1;
       x = Ntk + 1 + 3;
       
       % Абсолютні значення в напрямку польоту
       zvit{x, y} =...
           'Абсолютне відхилення в напрямку польоту, [м]';
       
       x = Ntk + 1 + 4;
       
       % Записуємо дані
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dabs3t{i, j}(it, ik);
          end
       end
       
       y = Ntk + 1 + 2;
       x = Ntk + 1 + 3;
       
       % Відносні значення в поперек польоту
       zvit{x, y} =...
           'Абсолютне відхилення в поперек польоту, [м]';
       
       x = Ntk + 1 + 4;
       
       % Записуємо дані
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dabs3k{i, j}(it, ik);
          end
       end
       
       % Запис у файл
       sheetname = strcat('Піксель (', string(i), ', ', string(j), ')');
       xlswrite(filename, zvit, sheetname);
   end
end

clear Dabs3k Dabs3t Dvid3k Dvid3t

% Аналізуємо дані
% Чистимо для подальших вивидень значень
zvit = {};

for it = 1:Ntk
   for ik = 1:Ntk
       theta = round((it - 1) * thetaE * (180.0 / pi), 5);
       phi = round((ik - 1) * phiE * (180.0 / pi), 5);
       
       y = 1;
       x = 1;
       
       % Відносні значення в напрямку польоту
       zvit{x, y} =...
           'Середнє ариф. відносне відхилення в напрямку польоту, [%]';

       x = 2;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaVid3t{it, ik}(:)) / (Nx * Ny);
       
       y = Ntk + 1 + 2;
       x = 1;
       
       % Відносні значення в поперек польоту
       zvit{x, y} =...
           'Середнє ариф. відносне відхилення в поперек польоту, [%]';
       
       x = 2;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaVid3k{it, ik}(:)) / (Nx * Ny);
       
       y = 1;
       x = Ntk + 1 + 3;
       
       % Абсолютні значення в напрямку польоту
       zvit{x, y} =...
           'Середнє ариф. абсолютне відхилення в напрямку польоту, [м]';
       
       x = Ntk + 1 + 4;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaAbs3t{it, ik}(:)) / (Nx * Ny);
       
       y = Ntk + 1 + 2;
       x = Ntk + 1 + 3;
       
       % Відносні значення в поперек польоту
       zvit{x, y} =...
           'Середнє ариф. абсолютне відхилення в поперек польоту, [м]';
       
       x = Ntk + 1 + 4;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaAbs3k{it, ik}(:)) / (Nx * Ny);
       
       % Загальне відносне відхилення
       y = 2 * (Ntk + 1 + 2) - 1;
       x = 1;
       
       % Відносні значення в напрямку польоту
       zvit{x, y} =...
           'Середнє геом. відносне відхилення, [%]';
       
       x = 2;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
%        zvit{it + x, ik + y} = ...
%            sqrt(sum(deltaVid3t{it, ik}(:)) / (Nx * Ny) *...
%            sum(deltaVid3k{it, ik}(:)) / (Nx * Ny));
       zvit{it + x, ik + y} = ...
           sqrt(sum(deltaVid3t{it, ik}(:)) *...
           sum(deltaVid3k{it, ik}(:))) / (Nx * Ny);
       
       y = 2 * (Ntk + 1 + 2) - 1;
       x = Ntk + 1 + 3;
       
       % Відносні значення в поперек польоту
       zvit{x, y} =...
           'Середнє геом. абсолютне відхилення, [м]';
       
       x = Ntk + 1 + 4;
       
       % Записуємо дані
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sqrt(sum(deltaAbs3t{it, ik}(:)) *...
           sum(deltaAbs3k{it, ik}(:))) / (Nx * Ny);
       
   end
end

clear deltaAbs3k deltaAbs3t deltaVid3k deltaVid3t phi theta phiE thetaE

% Запис у файл
xlswrite(filename, zvit, 'Аналіз');

clear Ntk i j x y ik it

% Виводимо сповіщення про збереження звіту
msgbox({'Розрахунок відхилень збережений у файлі:'...
    filename}, 'Звіт', 'help');

clear filename Phi Theta Nx Ny sheetname zvit
