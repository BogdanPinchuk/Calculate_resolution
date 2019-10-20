script % Projection_of_pixels
% ������ ������ F5
% ������ ��������� F9

%%%%%%%%
% ���� %
%%%%%%%%
% tic;
clear; clc;

% ʳ������ �������� �������� ��������
% � �������� ������� (�� ��������, Ox) TDI
pd = 8; % 32;
% ������� ������� (�� ����������, Oy)
qd = 32; % 4096;

% ����� ������

% ����� ��������� ��������
% ������� [���]
vd = 10; %8 * 10 ^ -6;
% ������ [���]
wd = 10; %8 * 10 ^ -6;

% ����� �������� ��������
% ������ [���]
Vd = 12; %8 * 10 ^ -6;
% �������� [���]
Wd = 12; %8 * 10 ^ -6;

% ������� ������ ���� ��������� � ����� �� (��� � �����)
% ������ [���]
dLp = 0; % 0 * 10 ^ -6;
% �������� [���]
dLq = 0; % 0 * 10 ^ -6;

% ������� ������� ��'������ [��]
f = 200; % 2260 * 10 ^ -3;

% ������ �����  [��]
h = 350; %668 * 10 ^ 3;

% ���� ��������
% ������� [����]
theta = 25 * (pi / 180.0);
% ����� [����]
phi = 25 * (pi / 180.0);
% �������� [����]
psi = 0 * (pi / 180.0);

% id - ������ ���� ������� �� ������������ �������
% id = 0 - �� �������� �� �������, � ��� �� �����;
% id = 1 - �� �������� �� �����, � ��� �� �������;
% id = 0;

%%%%%%%%%%%%%%
% ���������� %
%%%%%%%%%%%%%%

% theta_2 = atan(tan(phi) * cos(theta));
% phi_2 = atan(tan(theta) * cos(phi));
% alhpa = acos(cos(phi) * cos(phi_2));

% ��������� ����� ��������
Lp0 = dLp - 0.5 * Vd * (pd - 1);
Lq0 = dLq - 0.5 * Wd * (qd - 1);

clear dLp dLq

% ���������� ������ ��������� ������ ������
for i = 1:pd
    Lpi(i) = Lp0 + (i - 1) * Vd;
end

for i = 1:qd
   Lqi(i) = Lq0 + (i - 1) * Wd;
end

clear Vd Wd Lp0 Lq0 i

% ���������� ������� ��������� ������ ������
% �� �������� ��� ����� 0 (0, 0)
Wxi{1} = atan(Lpi ./ f);
Wyi{1} = atan(Lqi ./ f);

clear Lpi Lqi

% ����� �������� ������
wx = acot((f ./ (vd .* (cos(Wxi{1}) .^ 2))) - (vd / (4 * f)));
wy = acot((f ./ (wd .* (cos(Wyi{1}) .^ 2))) - (wd / (4 * f)));

% г����� ������� ������� ������ ������� �� ������
dwx = sign(Wxi{1}) .* acos(((vd / (2 * f)) .* (cos(Wxi{1}) .^ 2)...
    .* sin(wx) + cos(wx)));

dwy = sign(Wyi{1}) .* acos(((wd / (2 * f)) .* (cos(Wyi{1}) .^ 2)...
    .* sin(wy) + cos(wy)));

clear vd wd f

% �� �������� ��� ����� 1 (1, -1)
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);

% �� �������� ��� ����� 2 (1, 1)
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);

% �� �������� ��� ����� 3 (-1, 1)
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);

% �� �������� ��� ����� 1 (-1, -1)
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);

clear wx wy mx my dwx dwy

% � ��������� �������� ����� ����� ��� ����� ������
% ������� ����� ���� ��������
for k = 1:5
    for i = 1:pd
       for j = 1:qd
           theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
               - tan(Wyi{k}(1, j)) * sin(psi));
           phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
               + tan(Wxi{k}(1, i)) * sin(psi));
       end
    end
end

clear psi k i j Wxi Wyi

% ���������� ���� ������� ��������� � ��������� �� ���� ������
for k = 1:5
    % �������� �� �������, � ��� �� �����
    id = 0;
    [theta1{k}, phi1{k}] = naklon(theta, phi, theta0{k}, phi0{k}, id);
    
    % �������� �� �����, � ��� �� �������
    id = 1;
    [theta2{k}, phi2{k}] = naklon(theta, phi, theta0{k}, phi0{k}, id);
end

clear id k theta phi theta0 phi0

% ����������� ���� ���������� �����
for k = 1:5
    Mtheta1{k} = h .* tan(theta1{k});
    Mtheta2{k} = h .* tan(theta2{k});
    
    Mphi1{k} = h .* tan(phi1{k});
    Mphi2{k} = h .* tan(phi2{k});
end

clear h k theta1 phi1 theta2 phi2

% ������� ���������� ��������
figure(1);
clf;
hold on;
grid off;
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (������ - ����)');
% plot(0,0, 'k^');

% ��������� �������� � ����� ��� ���������� �������� �� ����
xmin = min(Mphi1{1}(:));
% ymin = min(Mtheta1{1}(:));
xmax = max(Mphi1{1}(:));
% ymax = max(Mtheta1{1}(:));
for k = 2:5
    xmin = min(min(Mphi1{k}(:)), xmin);
    xmax = max(max(Mphi1{k}(:)), xmax);
    xmin = min(min(Mtheta1{k}(:)), xmin);
    xmax = max(max(Mtheta1{k}(:)), xmax);
end

for k = 1:5
    xmin = min(min(Mphi2{k}(:)), xmin);
    xmax = max(max(Mphi2{k}(:)), xmax);
    xmin = min(min(Mtheta2{k}(:)), xmin);
    xmax = max(max(Mtheta2{k}(:)), xmax);
end

kl = 1.0;

xlim([kl * xmin kl * xmax]);
ylim([kl * xmin kl * xmax]);

% clear xmin ymin xmax ymax kl

%  ³��������� �������� ������
    for i = 1:pd
       for j = 1:qd
           plot(Mphi1{1}(i, j), Mtheta1{1}(i, j), 'r+');
           plot([Mphi1{2}(i, j), Mphi1{3}(i, j),...
               Mphi1{4}(i, j), Mphi1{5}(i, j), Mphi1{2}(i, j)],...
               [Mtheta1{2}(i, j), Mtheta1{3}(i, j), Mtheta1{4}(i, j),...
               Mtheta1{5}(i, j), Mtheta1{2}(i, j)], '-k');
       end
    end

figure(2);
clf;
hold on;
grid off;
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (���� - ������)');

% ��������� �������� � ����� ��� ���������� �������� �� ����
% xmin = min(Mphi2{1}(:));
% % ymin = min(Mtheta2{1}(:));
% xmax = max(Mphi2{1}(:));
% % ymax = max(Mtheta2{1}(:));
% for k = 2:5
%     xmin = min(min(Mphi2{k}(:)), xmin);
%     xmax = max(max(Mphi2{k}(:)), xmax);
%     xmin = min(min(Mtheta2{k}(:)), xmin);
%     xmax = max(max(Mtheta2{k}(:)), xmax);
% end
% 
% kl = 1.1;

xlim([kl * xmin kl * xmax]);
ylim([kl * xmin kl * xmax]);

clear xmin ymin xmax ymax kl

%  ³��������� �������� ������
    for i = 1:pd
       for j = 1:qd
           plot(Mphi2{1}(i, j), Mtheta2{1}(i, j), 'r+');
           plot([Mphi2{2}(i, j), Mphi2{3}(i, j),...
               Mphi2{4}(i, j), Mphi2{5}(i, j), Mphi2{2}(i, j)],...
               [Mtheta2{2}(i, j), Mtheta2{3}(i, j), Mtheta2{4}(i, j),...
               Mtheta2{5}(i, j), Mtheta2{2}(i, j)], '-k');
       end
    end

clear k i j pd qd
    
% ���������� ����� �� ������������� �������� �������� �������
% �������� �� �������, � ��� �� �����
Mtheta3{1} = centertoch(Mphi1{1}, Mtheta1{2}, Mtheta1{3},...
    Mphi1{2}, Mphi1{3}, 0);
Mphi3{1} = Mphi1{1};

Mtheta3{2} = Mtheta1{1};
Mphi3{2} = centertoch(Mtheta1{1}, Mtheta1{3}, Mtheta1{4},...
    Mphi1{3}, Mphi1{4}, 1);

Mtheta3{3} = centertoch(Mphi1{1}, Mtheta1{4}, Mtheta1{5},...
    Mphi1{4}, Mphi1{5}, 0);
Mphi3{3} = Mphi1{1};

Mtheta3{4} = Mtheta1{1};
Mphi3{4} = centertoch(Mtheta1{1}, Mtheta1{2}, Mtheta1{5},...
    Mphi1{2}, Mphi1{5}, 1);

% �������� �� �����, � ��� �� �������
Mtheta4{1} = centertoch(Mphi2{1}, Mtheta2{2}, Mtheta2{3},...
    Mphi2{2}, Mphi2{3}, 0);
Mphi4{1} = Mphi2{1};

Mtheta4{2} = Mtheta2{1};
Mphi4{2} = centertoch(Mtheta2{1}, Mtheta2{3}, Mtheta2{4},...
    Mphi2{3}, Mphi2{4}, 1);

Mtheta4{3} = centertoch(Mphi2{1}, Mtheta2{4}, Mtheta2{5},...
    Mphi2{4}, Mphi2{5}, 0);
Mphi4{3} = Mphi2{1};

Mtheta4{4} = Mtheta2{1};
Mphi4{4} = centertoch(Mtheta2{1}, Mtheta2{2}, Mtheta2{5},...
    Mphi2{2}, Mphi2{5}, 1);

% ����������� ���������� �������� �������� ������ �������
% �������� �� �������, � ��� �� �����
Ltheta1 = centertoch(0, Mtheta3{1}, Mtheta3{3},...
    Mphi3{1}, Mphi3{3}, 2);

Lphi1 = centertoch(0, Mtheta3{2}, Mtheta3{4},...
    Mphi3{2}, Mphi3{4}, 2);

% �������� �� �����, � ��� �� �������
Ltheta2 = centertoch(0, Mtheta4{1}, Mtheta4{3},...
    Mphi4{1}, Mphi4{3}, 2);

Lphi2 = centertoch(0, Mtheta4{2}, Mtheta4{4},...
    Mphi4{2}, Mphi4{4}, 2);

clear Mphi1 Mphi2 Mtheta1 Mtheta2 Mphi3 Mphi4 Mtheta3 Mtheta4

% % ���������� ���� �������� ������
% % ������� ���������
% % �������� �� �������, � ��� �� �����
% Dvz1 = centertoch(0, Mtheta1{2}, Mtheta1{4},...
%     Mphi1{2}, Mphi1{4}, 2);
% Dpo1 = centertoch(0, Mtheta1{3}, Mtheta1{5},...
%     Mphi1{3}, Mphi1{5}, 2);
% 
% % �������� �� �����, � ��� �� �������
% Dvz2 = centertoch(0, Mtheta2{2}, Mtheta2{4},...
%     Mphi2{2}, Mphi2{4}, 2);
% Dpo2 = centertoch(0, Mtheta2{3}, Mtheta2{5},...
%     Mphi2{3}, Mphi2{5}, 2);
% 
% % �������� ����������� ��� ���������� ���� �� ����������
% % �������� �� �������, � ��� �� �����
% Lt31_1 = Mtheta1{4} - Mtheta1{2};
% Lt42_1 = Mtheta1{5} - Mtheta1{3};
% Lp31_1 = Mphi1{4} - Mphi1{2};
% Lp42_1 = Mphi1{5} - Mphi1{3};
% 
% % �������� �� �����, � ��� �� �������
% Lt31_2 = Mtheta2{4} - Mtheta2{2};
% Lt42_2 = Mtheta2{5} - Mtheta2{3};
% Lp31_2 = Mphi2{4} - Mphi2{2};
% Lp42_2 = Mphi2{5} - Mphi2{3};
 

% ��� ���������� [�]/[��]
% t = toc; % / 60.0;
% clear t
clc;

% for i = 1:Nx
%     if (tan(phi1{1}(i, 1)) - tan(phi1{1}(i, Ny))) ~= 0
%         PhiR1(i) = atan((tan(theta1{1}(i, 1)) - tan(theta1{1}(i, Ny))) /...
%             (tan(phi1{1}(i, 1)) - tan(phi1{1}(i, Ny))));
%         % � ��������
%         PhiR1(i) = 90 - PhiR1(i) * (180 / pi);
%     else
%         % � ��������
%         PhiR1(i) = 90;
%     end
%     
%     if (tan(phi2{1}(i, 1)) - tan(phi2{1}(i, Ny))) ~= 0
%         PhiR2(i) = atan((tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) /...
%             (tan(phi2{1}(i, 1)) - tan(phi2{1}(i, Ny))));
%         % � ��������
%         PhiR2(i) = 90 + PhiR2(i) * (180 / pi);
%     else
%         % � ��������
%         PhiR2(i) = 90;
%     end
% end

% for j = 1:Ny
%     if (tan(phi1{1}(1, j)) - tan(phi1{1}(Nx, j))) ~= 0
%         PhiC1(j) = atan((tan(theta1{1}(1, j)) - tan(theta1{1}(Nx, j))) /...
%             (tan(phi1{1}(1, j)) - tan(phi1{1}(Nx, j))));
%         % � ��������
%         PhiC1(j) = PhiC1(j) * (180 / pi);
%     else
%         % � ��������
%         PhiC1(j) = 90;
%     end
%     
%     if (tan(phi2{1}(1, j)) - tan(phi2{1}(Nx, j))) ~= 0
%         PhiC2(j) = atan((tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) /...
%             (tan(phi2{1}(1, j)) - tan(phi2{1}(Nx, j))));
%         % � ��������
%         PhiC2(j) = PhiC2(j) * (180 / pi);
%     else
%         % � ��������
%         PhiC2(j) = 90;
%     end
% end
