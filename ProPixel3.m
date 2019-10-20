script % Projection_of_pixels_2
% ������ ������ F5
% ������ ��������� F9

%%%%%%%%
% ���� %
%%%%%%%%
clear; clc;

% ʳ������ �������� �������� ��������
% � �������� ������� (�� ��������, Ox) TDI
pd = 32 + 1;
% ������� ������� (�� ����������, Oy)
qd = 4096 + 1;

% ����� ������

% ����� ��������� ��������
% ������� [���]
vd = 8 * 10 ^ -6;
% ������ [���]
wd = 8 * 10 ^ -6;

% ����� �������� ��������
% ������ [���]
Vd = 8 * 10 ^ -6;
% �������� [���]
Wd = 8 * 10 ^ -6;

% ������� ������� ��'������ [��]
f = 2260 * 10 ^ -3;

% ������� ������ ���� ��������� � ����� �� (��� � �����)
% ������ [����]
dAp = 0 * (pi / 180.0);
% �������� [����]
dAq = 0 * (pi / 180.0);

% ������ �����  [��]
h = 668 * 10 ^ 3;

% ���� ��������
% ������� [����]
theta = 0 * (pi / 180.0);
% ����� [����]
phi = 0 * (pi / 180.0);
% �������� [����]
psi = 0 * (pi / 180.0);

%%%%%%%%%%%%%%%%%%%%%%%%
% ��������� ��������� %
%%%%%%%%%%%%%%%%%%%%%%%%

% ʳ������ ������� ������ ��� �����������
% ������ �������� �������
Nx0 = 7;

% ������� �������� �������
Ny0 = 11;

% ³���������� �����
% ���� id = 0 - ���������� ��
% ���� id = 1 - ���������� ����� ������
idx = 1;
idy = 1;
% ���������� �������������
id = 1;

% ������� �� �������� � % (��������)
proc = 10;

%%%%%%%%%%%%%%
% ���������� %
%%%%%%%%%%%%%%

% �������� ������� ����� �� ��������� ����������
Nx = numbertoch(Nx0, pd, idx);
Ny = numbertoch(Ny0, qd, idy);

% ������� ������ ���� ��������� � ����� �� (��� � �����)
% ������ [���]
dLp = f * tan(dAp);
% �������� [���]
dLq = f * tan(dAq);

clear dAp dAq

% ��������� ����� ��������
Lp0 = dLp - 0.5 * Vd * (pd - 1);
Lq0 = dLq - 0.5 * Wd * (qd - 1);

clear dLp dLq

% ��������� ��������� ������������� ������
Mx = masivtoch(Nx, pd, idx);
My = masivtoch(Ny, qd, idy);

% �������� ������� ����� �� ������� ����������
Nx = length(Mx);
Ny = length(My);

clear pd qd

% ���������� ������ ��������� ������ ������
for i = 1:Nx
    Lpi(i) = Lp0 + (Mx(i) - 1) * Vd;
end

for i = 1:Ny
    Lqi(i) = Lq0 + (My(i) - 1) * Wd;
end

clear i Lp0 Lq0

% ���������� ������� ��������� ������ ������
% �� �������� ��� ����� "0" (0, 0)
Wxi{1} = atan(Lpi ./ f);
Wyi{1} = atan(Lqi ./ f);

% ���������� ������ ����� ��� ������ (�������)
WXi{1} = Wxi{1};
WYi{1} = Wyi{1};

clear Lpi Lqi

% ������� ����� ������� ������� ������ ��������������� ���� ���
% ����������� ��������, � ���������� �������� �������� �����������
% ��� ������� ������, �� ��������� ������������� 2 ������; ����� �������
% �������, �� 2-� ����� � ����������� ������ ������ ���������������������
% ��� ��������� ����������� ������ (�������� ����������)

% ����� �������� ������� ������� (��) ������
wx = anglesize(Wxi{1}, vd, f);
wy = anglesize(Wyi{1}, wd, f);

% ����� �������� ������ (�������)
Wx = anglesize(Wxi{1}, Vd, f);
Wy = anglesize(Wyi{1}, Wd, f);

% г����� ������� ������� �� ������ ������� �� ������
dwx = diffangle(Wxi{1}, wx, vd, f);
dwy = diffangle(Wyi{1}, wy, wd, f);

% г����� ������� ������� ������ (�������) ������� �� ������
dWx = diffangle(Wxi{1}, Wx, Vd, f);
dWy = diffangle(Wyi{1}, Wy, Wd, f);

clear vd wd Vd Wd f

% �� �������� ��� ����� 1 (1, -1)
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{2} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 2 (1, 1)
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{3} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 3 (-1, 1)
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{4} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 1 (-1, -1)
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{5} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{5} = koordtoch(Wyi{1}, Wy, my, dWy);

clear mx my wx wy dwx dwy Wx Wy dWx dWy

% � ��������� �������� ����� ����� ��� ����� ������
% ������� ����� ���� ��������
for k = 1:5
    for i = 1:Nx
       for j = 1:Ny
           % ��� �� ������
           theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
               - tan(Wyi{k}(1, j)) * sin(psi));
           phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
               + tan(Wxi{k}(1, i)) * sin(psi));
           
           % ��� �񳺿 ������� (��) ������
           Theta0{k}(i, j) = atan(tan(WXi{k}(1, i)) * cos(psi)...
               - tan(WYi{k}(1, j)) * sin(psi));
           Phi0{k}(i, j) = atan(tan(WYi{k}(1, j)) * cos(psi)...
               + tan(WXi{k}(1, i)) * sin(psi));
       end
    end
end

clear psi k i j Wxi Wyi WXi WYi

% ���������� ���� ������� ��������� � ��������� �� ���� ������
for k = 1:5
    % �������� �� �������, � ��� �� �����
    [theta1{k}, phi1{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 0);
    [Theta1{k}, Phi1{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 0);
    
    % �������� �� �����, � ��� �� �������
    [theta2{k}, phi2{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 1);
    [Theta2{k}, Phi2{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 1);
end

clear k theta phi theta0 phi0 Theta0 Phi0

% ���������� ���� ������ ������� � ��������
% ��� ����� (������� �������)
for i = 1:Nx
    if (tan(theta1{1}(i, 1)) - tan(theta1{1}(i, Ny))) ~= 0
        PhiR1(i) = acot((tan(theta1{1}(i, 1)) - tan(theta1{1}(i, Ny))) \...
            (tan(phi1{1}(i, 1)) - tan(phi1{1}(i, Ny))));
        % � ��������
        PhiR1(i) = 90 - PhiR1(i) * (180 / pi);
    else
        % � ��������
        PhiR1(i) = 90;
    end
    
    if (tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) ~= 0
        PhiR2(i) = acot((tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) \...
            (tan(phi2{1}(i, 1)) - tan(phi2{1}(i, Ny))));
        % � ��������
        PhiR2(i) = 90 - PhiR2(i) * (180 / pi);
    else
        % � ��������
        PhiR2(i) = 90;
    end
end

clear i

% ��� �������� (� �������� �������)
for j = 1:Ny
    if (tan(theta1{1}(1, j)) - tan(theta1{1}(Nx, j))) ~= 0
        PhiC1(j) = atan((tan(theta1{1}(1, j)) - tan(theta1{1}(Nx, j))) \...
            (tan(phi1{1}(1, j)) - tan(phi1{1}(Nx, j))));
        % � ��������
        PhiC1(j) = PhiC1(j) * (180 / pi);
    else
        % � ��������
        PhiC1(j) = 90;
    end
    
    if (tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) ~= 0
        PhiC2(j) = atan((tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) \...
            (tan(phi2{1}(1, j)) - tan(phi2{1}(Nx, j))));
        % � ��������
        PhiC2(j) = PhiC2(j) * (180 / pi);
    else
        % � ��������
        PhiC2(j) = 90;
    end
end

clear j

% ����������� ���� ���������� �����
for k = 1:5
    % ��� �� ������
    Mtheta1{k} = h .* tan(theta1{k});
    Mtheta2{k} = h .* tan(theta2{k});
    
    Mphi1{k} = h .* tan(phi1{k});
    Mphi2{k} = h .* tan(phi2{k});
    
    % ��� �� ������
    MTheta1{k} = h .* tan(Theta1{k});
    MTheta2{k} = h .* tan(Theta2{k});
    
    MPhi1{k} = h .* tan(Phi1{k});
    MPhi2{k} = h .* tan(Phi2{k});
end

clear h k Theta1 Phi1 Theta2 Phi2 theta1 phi1 theta2 phi2 


% ³���������� ����� ���������� ������, � �� ������ ������
if (id == 1) && (idx == 1) && (idy == 1) && (Nx0 > 1) && (Ny0 > 1)
    % �������� �ٲ�������
    % ���������� ��������� ����� � ����������� �����
    Mx = sysidtoch(Mx);
    My = sysidtoch(My);

    % ��������� ����� ��� ��� ������ (center weight - CW)
    [CWx1, CWy1] = centerweight(Mx, My, MTheta1, MPhi1, Nx, Ny);
    [CWx2, CWy2] = centerweight(Mx, My, MTheta2, MPhi2, Nx, Ny);

    % ��������� ������� ��� ��������� ��������� ������� �����������
    [MinT1, MinP1] = minvidmizpix(Mx, My, CWx1, CWy1, Nx, Ny, proc);
    [MinT2, MinP2] = minvidmizpix(Mx, My, CWx2, CWy2, Nx, Ny, proc);
    
    % ���������� �� ���������� ������������ ����� ������
    [Mtheta1, Mphi1] = everyscale(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    % ��� ���������� ���� �� ������� 1 ��������������, ��� �������� ���
    % ����� ��� �������, ������� ��������� ��������� ��� ������ ���
    [CWx1, CWy1] = centerweight2(MTheta1, MPhi1, MinT1, MinP1,...
        Nx, Ny, CWx1, CWy1);
    [CWx2, CWy2] = centerweight2(MTheta2, MPhi2, MinT2, MinP2,...
        Nx, Ny, CWx2, CWy2);
    
    % ��������� �������������� ������� ���, ��� ���������� �������������
    % �� ����� �� ����
    % ��������� ����������� ���� �������������
    Min1 = min(MinT1, MinP1);
    MinT1 = MinT1 / Min1;
    MinP1 = MinP1 / Min1;
    Min2 = min(MinT2, MinP2);
    MinT2 = MinT2 / Min2;
    MinP2 = MinP2 / Min2;
    
    clear Min1 Min2
    
    % ���������� �� ���������� ������������ ����� ������
    [Mtheta1, Mphi1] = everyscale2(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale2(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    clear MinT1 MinP1 MinT2 MinP2 %CWx1 CWy1 CWx2 CWy2
end

clear idx idy Nx0 Ny0 Mx My

% ������� ���������� ��������
figure(1);
clf;
hold on;
grid off;   % �������� ���� "on", ��������� "off"
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (������ - ����)');
% plot(0,0, 'k^');

% ��������� �������� � ����� ��� ���������� �������� �� ����
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

% ���������� ���������� ��
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

%  ³��������� �������� ������
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

% �������� �� ������ ������ ������ �� ���� ��
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
grid off;   % �������� ���� "on", ��������� "off"
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (���� - ������)');

% ylim([ymin ymax]);
% xlim([xmin xmax]);

ylim(xcenter - 0.5 .* [len -len]);
xlim(ycenter - 0.5 .* [len -len]);

clear xcenter ycenter len

%  ³��������� �������� ������
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

% �������� �� ������ ������ ������ �� ���� ��
% plot([Mphi2{1}(Nx, 1), Mphi2{1}(Nx, Ny),],...
%     [Mtheta2{1}(Nx, 1), Mtheta2{1}(Nx, Ny)], '-g');
% plot([Mphi2{1}(1, 1), Mphi2{1}(Nx, 1),],...
%     [Mtheta2{1}(1, 1), Mtheta2{1}(Nx, 1)], '-g');
% 
% plot([Mphi2{1}(1, 1), Mphi2{1}(1, Ny),],...
%     [Mtheta2{1}(1, 1), Mtheta2{1}(1, Ny)], '-g');

clear k i j pd qd id CWx1 CWy1 CWx2 CWy2 Mphi1 Mphi2 Mtheta1 Mtheta2

% ���������� ����� �� ������������� �������� �������� �������
% ����� �������� �� �������, � ��� �� �����
[Ltheta1, Lphi1] = razreshtoch(MTheta1, MPhi1, Nx, Ny);
[Ltheta2, Lphi2] = razreshtoch(MTheta2, MPhi2, Nx, Ny);

clear MTheta1 MPhi1 MTheta2 MPhi2 Nx Ny

% ���������� ������� [�� ^ -1]
NuTheta1 = 1000 ./ Ltheta1;
NuLphi1 = 1000 ./ Lphi1;
NuTheta2 = 1000 ./ Ltheta2;
NuLphi2 = 1000 ./ Lphi2;

% figure(1);
% hold on;
% %  ³��������� �������� ������
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
% %  ³��������� �������� ������
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
