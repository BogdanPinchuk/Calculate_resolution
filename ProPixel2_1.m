script % Projection_of_pixels_2
% ������ ������ F5
% ������ ��������� F9

%%%%%%%%
% ���� %
%%%%%%%%
clear; clc;

% ʳ������ �������� �������� ��������
% � �������� ������� (�� ��������, Ox) TDI
pd = 33; % 32;
% ������� ������� (�� ����������, Oy)
qd = 163; % 4096;

% ����� ������

% ����� ��������� ��������
% ������� [���]
vd = 13 * 10 ^ -6; %8 * 10 ^ -6;
% ������ [���]
wd = 13 * 10 ^ -6; %8 * 10 ^ -6;

% ����� �������� ��������
% ������ [���]
Vd = 13 * 10 ^ -6; %8 * 10 ^ -6;
% �������� [���]
Wd = 13 * 10 ^ -6; %8 * 10 ^ -6;

% ������� ������� ��'������ [��]
f = 9650 * 10 ^ -3; % 2260 * 10 ^ -3;

% ������� ������ ���� ��������� � ����� �� (��� � �����)
% ������ [����]
dAp = 0 * (pi / 180.0);
% �������� [����]
dAq = 0 * (pi / 180.0);

% ������ �����  [��]
h = 668 * 10 ^ 3;

% ���� ��������
% ������� [����]
theta = 25 * (pi / 180.0);
% ����� [����]
phi = 25 * (pi / 180.0);
% �������� [����]
psi = 0 * (pi / 180.0);

% ������ (-90�...+90�) �� ��� ����������� �� [����]
gama = 50 * (pi / 180.0);

% ������ ����, ����� ������� �� ����� ������� ���� �� �� ������ [��]
Rzmin = 6356.777 * 10 ^ 3;
Rzmax = 6378.160 * 10 ^ 3;

%%%%%%%%%%%%%%%%%%%%%%%%
% ��������� ��������� %
%%%%%%%%%%%%%%%%%%%%%%%%

% ʳ������ ������� ������ ��� �����������
% ������ �������� �������
Nx0 = 5;

% ������� �������� �������
Ny0 = 5;

% ³���������� �����
% ���� id = 0 - ���������� ��
% ���� id = 1 - ���������� ����� ������
idx = 1;
idy = 1;
% ���������� ������������� (0 - �, 1 - ���)
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
Lp0 = 0.5 * Vd * (1 - pd) - dLp;
Lq0 = 0.5 * Wd * (1 - qd) - dLq;

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

% ���������� ������ ����� ��� ����. ����. ��
% wxi{1} = Wxi{1};
% wyi{1} = Wyi{1};

clear Lpi Lqi

% ������� ����� ������� ������� ������ ��������������� ���� ���
% ����������� ��������, � ���������� �������� �������� �����������
% ��� ������� ������, �� ��������� ������������� 2 ������; ����� �������
% ���������, �� 2-� ����� � ����������� ������ ������ ���������������������
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

clear vd wd

% ���������� ������� ������� �� (Wxi, Wyi)
% ���������� �񳺿 ������� �� (WXi, Wyi)
% ���������� ��� �������� �������� (wxi, wyi)
% �� �������� ��� ����� 1 (1, -1) �� (1, 0) 5
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
wxi{1} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{2} = koordtoch(Wyi{1}, Wy, my, dWy);
my = 0;
wyi{1} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 2 (1, 1) �� (0, 1) 6
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
mx = 0;
wxi{2} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{3} = koordtoch(Wyi{1}, Wy, my, dWy);
wyi{2} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 3 (-1, 1) �� (-1, 0) 7
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
wxi{3} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{4} = koordtoch(Wyi{1}, Wy, my, dWy);
my = 0;
wyi{3} = koordtoch(Wyi{1}, Wy, my, dWy);

% �� �������� ��� ����� 1 (-1, -1) �� (0, -1) 8
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
WXi{5} = koordtoch(Wxi{1}, Wx, mx, dWx);
mx = 0;
wxi{4} = koordtoch(Wxi{1}, Wx, mx, dWx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);
WYi{5} = koordtoch(Wyi{1}, Wy, my, dWy);
wyi{4} = koordtoch(Wyi{1}, Wy, my, dWy);

clear mx my wx wy dwx dwy dWx dWy %Wx Wy

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
           
           % ��� �������� ��������
           if k < 5
               Rtheta0{k}(i, j) = atan(tan(wxi{k}(1, i)) * cos(psi)...
                   - tan(wyi{k}(1, j)) * sin(psi));
               Rphi0{k}(i, j) = atan(tan(wyi{k}(1, j)) * cos(psi)...
                   + tan(wxi{k}(1, i)) * sin(psi));
           end
       end
    end
end

clear psi k i j Wxi Wyi WXi WYi wxi wyi

% ���������� ���� ������� ��������� � ��������� �� ���� ������
for k = 1:5
    % �������� �� �������, � ��� �� �����
    [theta1{k}, phi1{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 0);
    [Theta1{k}, Phi1{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 0);

    % �������� �� �����, � ��� �� �������
    [theta2{k}, phi2{k}] = naklon(theta, phi, theta0{k}, phi0{k}, 1);
    [Theta2{k}, Phi2{k}] = naklon(theta, phi, Theta0{k}, Phi0{k}, 1);
    
    if k < 5
        % �������� �� �������, � ��� �� �����
        [Rtheta1{k}, Rphi1{k}] = naklon(theta, phi, Rtheta0{k}, Rphi0{k}, 0);
        % �������� �� �����, � ��� �� �������
        [Rtheta2{k}, Rphi2{k}] = naklon(theta, phi, Rtheta0{k}, Rphi0{k}, 1);
    end
end

clear k theta phi theta0 phi0 Theta0 Phi0 Rtheta0 Rphi0

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
    
    % ��� �� �������� ��������
    if k < 5
        MRtheta1{k} = h .* tan(Rtheta1{k});
        MRtheta2{k} = h .* tan(Rtheta2{k});
        
        MRphi1{k} = h .* tan(Rphi1{k});
        MRphi2{k} = h .* tan(Rphi2{k});
    end
end

clear k Theta1 Phi1 Theta2 Phi2 theta1 phi1 theta2 phi2

% ���������� "������ �������" (������������� �� ������� ������ ������)
[LzahT1, LzahP1, CzahT1, CzahP1] = poloszah(Mtheta1, Mphi1, Nx, Ny);
[LzahT2, LzahP2, CzahT2, CzahP2] = poloszah(Mtheta2, Mphi2, Nx, Ny);

% ³���������� ����� ���������� ������, � �� ������ ������
if (id == 1) && (idx == 1) && (idy == 1) && (Nx0 > 1) && (Ny0 > 1)
    % �������� �ٲ�������
    % ���������� ��������� ����� � ����������� �����
    Kx = sysidtoch(Mx);
    Ky = sysidtoch(My);
    
    clear Mx My
    
    % ��������� ����� ��� ��� ������ (center weight - CW)
    [CWx1, CWy1] = centerweight(Kx, Ky, MTheta1, MPhi1, Nx, Ny);
    [CWx2, CWy2] = centerweight(Kx, Ky, MTheta2, MPhi2, Nx, Ny);

    % ��������� ������� ��� ��������� ��������� ������� �����������
    [MinT1, MinP1] = minvidmizpix(Kx, Ky, CWx1, CWy1, Nx, Ny, proc);
    [MinT2, MinP2] = minvidmizpix(Kx, Ky, CWx2, CWy2, Nx, Ny, proc);
    
    % ��������� ������������ ����. ��������, ��� ���� ��� ����� �������
    % �������������� �������� ������ � � ��������� ���� �� ������ ������
    if ~isnan(MinT1) && ~isnan(MinP1)
        MaxSc1 = max(MinT1, MinP1);
    else
        MaxSc1 = 1.0;
    end
        
    if ~isnan(MinT2) && ~isnan(MinP2)
        MaxSc2 = max(MinT2, MinP2);
    else
        MaxSc2 = 1.0;
    end
    
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
    
    % ����������� ��� ��� ����������� ����������
    minT1 = min(Mtheta1{1}(:));
    maxT1 = max(Mtheta1{1}(:));
    minP1 = min(Mphi1{1}(:));
    maxP1 = max(Mphi1{1}(:));
    
    minT2 = min(Mtheta2{1}(:));
    maxT2 = max(Mtheta2{1}(:));
    minP2 = min(Mphi2{1}(:));
    maxP2 = max(Mphi2{1}(:));
   
    for k = 2:5
        minT1 = min(minT1, min(Mtheta1{k}(:)));
        maxT1 = max(maxT1, max(Mtheta1{k}(:)));
        minP1 = min(minP1, min(Mphi1{k}(:)));
        maxP1 = max(maxP1, max(Mphi1{k}(:)));

        minT2 = min(minT2, min(Mtheta2{k}(:)));
        maxT2 = max(maxT2, max(Mtheta2{k}(:)));
        minP2 = min(minP2, min(Mphi2{k}(:)));
        maxP2 = max(maxP2, max(Mphi2{k}(:)));
    end
    
    % ��������� ������
    centT1 = minT1 + 0.5 * (maxT1 - minT1);
    centP1 = minP1 + 0.5 * (maxP1 - minP1);
    centT2 = minT2 + 0.5 * (maxT2 - minT2);
    centP2 = minP2 + 0.5 * (maxP2 - minP2);
    
    % �������� � ����������
    for k = 1:5
        Mtheta1{k} = (Mtheta1{k} - centT1) / MaxSc1;
        Mphi1{k} = (Mphi1{k} - centP1) / MaxSc1;
        Mtheta2{k} = (Mtheta2{k} - centT2) / MaxSc2;
        Mphi2{k} = (Mphi2{k} - centP2) / MaxSc2;
    end
    
    clear MinT1 MinP1 MinT2 MinP2... %CWx1 CWy1 CWx2 CWy2
    minT1 maxT1 minP1 maxP1 minT2 maxT2 minP2 maxP2...
    centT1 centP1 centT2 centP2 MaxSc1 MaxSc2
end

clear idx idy Nx0 Ny0 Kx Ky Mx My

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
       
       % �������� ���������� ���������� ����� ��� ��
       if id == 0
           for k = 1:4
               plot(MRphi1{k}(i, j), MRtheta1{k}(i, j), 'b+');
           end
       end
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
       
       % �������� ���������� ���������� ����� ��� ��
       if id == 0
           for k = 1:4
               plot(MRphi2{k}(i, j), MRtheta2{k}(i, j), 'b+');
           end
       end
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

clear k i j pd qd id CWx1 CWy1 CWx2 CWy2 %Mphi1 Mphi2 Mtheta1 Mtheta2

% ���������� ����� �� ������������� �������� �������� ������� (�� ������)
% 1-�� ������ (��� ���������� �������� � ������ �� ����)
% [Ltheta1, Lphi1] = razreshtoch(MTheta1, MPhi1, Nx, Ny);
% [Ltheta2, Lphi2] = razreshtoch(MTheta2, MPhi2, Nx, Ny);

% 2-�� ������ ���������� ������ �� �������� �������
% �������� �� �������, � ��� �� �����
Ltheta1 = vipstanmizpix(MRtheta1{1}, MRtheta1{3}, MRphi1{1}, MRphi1{3});
Lphi1 = vipstanmizpix(MRtheta1{2}, MRtheta1{4}, MRphi1{2}, MRphi1{4});

% �������� �� �����, � ��� �� �������
Ltheta2 = vipstanmizpix(MRtheta2{1}, MRtheta2{3}, MRphi2{1}, MRphi2{3});
Lphi2 = vipstanmizpix(MRtheta2{2}, MRtheta2{4}, MRphi2{2}, MRphi2{4});

clear MTheta1 MPhi1 MTheta2 MPhi2 MRtheta1 MRphi1 MRtheta2 MRphi2

% ���������� ��������� ������� [�� ^ -1]
NuTheta1 = (10 ^ -3) ./ Ltheta1;
NuPhi1 = (10 ^ -3) ./ Lphi1;
NuTheta2 = (10 ^ -3) ./ Ltheta2;
NuPhi2 = (10 ^ -3) ./ Lphi2;

% ��������� ������� 

%clear NuTheta1 NuTheta2 NuLphi1 NuLphi2 PhiC1 PhiC2 PhiR1 PhiR2...
    %Lphi1 Lphi2 Ltheta1 Ltheta2

% ���������� �������� �������� ���������� �������� ����
% ����������� ������� �� ������ �� ����� �� ������� ����
Rt = radzem(Rzmin, Rzmax, gama);

% ����������� ����� �������� ����
Rk = radkrivzem(Rzmin, Rzmax, gama);

clear gama

% ����������� ������� �� ������
% ������ ���� 45�, �� ��� ������� �� ������ � �� ��������� = h
dh = Rt - radzem(Rzmin, Rzmax, 45 * (pi / 180));

clear Rt Rzmin Rzmax

% ��������� ������
h = h + dh;

clear dh

% ����������� ���� ��������� �� ������
for i = 1:Nx
    for j = 1:Ny
        for k = 1:4
            alpha1{k}(i, j) = anglealpha(Rtheta1{k}(i, j), Rphi1{k}(i, j));
            alpha2{k}(i, j) = anglealpha(Rtheta2{k}(i, j), Rphi2{k}(i, j));
        end
    end
end

clear i j k Rtheta1 Rphi1 Rtheta2 Rphi2

% ��������� ������� �� �� �� ����� ����� ��� ����� �� ������
for k = 1:4
    dovzin1{k} = dovzdotoch(h, Rk, alpha1{k});
    dovzin2{k} = dovzdotoch(h, Rk, alpha2{k});
end

clear alpha1 alpha2 h k

% ��������� �������� �������� ���������� ��������, ��� ����������� ��
% �����
[LRatheta1, LRaphi1] = vidsmiztoch(dovzin1, Wx, Wy, Nx, Ny);
[LRatheta2, LRaphi2] = vidsmiztoch(dovzin2, Wx, Wy, Nx, Ny);

clear dovzin1 dovzin2

% ���������� ������� ������� ������
[WX, WY] = angleRZ(Wx, Wy, Vd, Wd, f, Nx, Ny);

% ���������� ����� ������� [���� ^ -1]
NuWx = (10 ^ -3) ./ WX;
NuWy = (10 ^ -3) ./ WY;

clear Vd Wd f Nx Ny Wx Wy WX WY

% % ���������� �������� �������� ���������� ��������
% LRtheta1 = realrozdzdat(LRatheta1, Rk);
% LRphi1 = realrozdzdat(LRaphi1, Rk);
% LRtheta2 = realrozdzdat(LRatheta2, Rk);
% LRphi2 = realrozdzdat(LRaphi2, Rk);

clear Rk

% г����� �� ������������� � �������� �� � ����������� ��������
% dLRtheta1 = 100 .* (LRatheta1 - Ltheta1) ./ Ltheta1;
% dLRphi1 = 100 .* (LRaphi1 - Lphi1) ./ Lphi1;
% dLRtheta2 = 100 .* (LRatheta2 - Ltheta2) ./ Ltheta2;
% dLRphi2 = 100 .* (LRaphi2 - Lphi2) ./ Lphi2;

% dlLRtheta1 = LRatheta1 - Ltheta1;
% dlLRphi1 = LRaphi1 - Lphi1;
% dlLRtheta2 = LRatheta2 - Ltheta2;
% dlLRphi2 = LRaphi2 - Lphi2;

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

% ���� ���������� ������� ����������

% roundn([min(dLRtheta1(:)) max(dLRtheta1(:)); min(dLRphi1(:)) max(dLRphi1(:))], -1)
% roundn([min(dLRtheta2(:)) max(dLRtheta2(:)); min(dLRphi2(:)) max(dLRphi2(:))], -1)

% roundn([min(dlLRtheta1(:)) max(dlLRtheta1(:)); min(dlLRphi1(:)) max(dlLRphi1(:))], -2)
% roundn([min(dlLRtheta2(:)) max(dlLRtheta2(:)); min(dlLRphi2(:)) max(dlLRphi2(:))], -2)

clear ans
