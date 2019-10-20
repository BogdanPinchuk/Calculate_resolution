script % Projection_of_pixels_4
% ������ ������ F5
% ������ ��������� F9

%%%%%%%%
% ���� %
%%%%%%%%
clear; clc;

% ʳ������ �������� �������� ��������
% � �������� ������� (�� ��������, Ox) TDI
pd = 41;
% ������� ������� (�� ����������, Oy)
qd = 40001;

% ����� ������

% ����� ��������� ��������
% ������� [���]
vd = 13 * 10 ^ -6;
% ������ [���]
wd = 13 * 10 ^ -6;

% ����� �������� ��������
% ������ [���]
Vd = 13 * 10 ^ -6;
% �������� [���]
Wd = 13 * 10 ^ -6;

% ������� ������� ��'������ [��]
f = 9650 * 10 ^ -3;

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

% ������ (-90�...+90�) �� ��� ����������� �� [����]
gama = 50 * (pi / 180.0);

% ������ ����, ����� ������� �� ����� ������� ���� �� �� ������ [��]
Rz = 6371.032 * 10 ^ 3;
Rzmin = 6356.777 * 10 ^ 3;
Rzmax = 6378.160 * 10 ^ 3;

%%%%%%%%%%%%%%%%%%%%%%%%
% ��������� ��������� %
%%%%%%%%%%%%%%%%%%%%%%%%

% ʳ������ ������� ������ ��� �����������
% ������ �������� �������
Nx0 = 7;

% ������� �������� �������
Ny0 = 7;

% ³���������� �����
% ���� id = 0 - ���������� ��
% ���� id = 1 - ���������� ����� ������
idx = 1;
idy = 1;
% ���������� ������������� (0 - �, 1 - ���)
id = 1;

% ������� �� �������� � % (��������)
proc = 10;

% ������� ���������� 1 (false) ��� 
% ����� ������ (true)
Npix = true;

% ���������� ��� ��� true ������� �� ����� false
CR = true;

% ���������� ������� true, ��� ����� false
HV = true;

% ����� ������� ��� �����
numHV = 20000;
 
 % ��������� ��������� ������ (���� �� �����, ����� �� �����)
pixNx = 20;
pixNy = 20000;
 
%%%%%%%%%%%%%%
% ���������� %
%%%%%%%%%%%%%%

% �������� ������� ����� �� ��������� ����������
if Npix % ���������� 1 ������ ��� ���
    if CR % ���������� ������ ��� ������ �������
        Nx = numbertoch(Nx0, pd, idx);
        Ny = numbertoch(Ny0, qd, idy);
    else
        if HV % ���������� �������
            Nx = numbertoch(Nx0, pd, idx);
            Ny = 1;
        else % ���������� �����
            Nx = 1;
            Ny = numbertoch(Ny0, qd, idy);
        end
    end
else
    Nx = 1;
    Ny = 1;
end

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
if Npix % ���������� 1 ������ ��� ���
    if CR % ���������� ������ ��� ������ �������
        Mx = masivtoch(Nx, pd, idx);
        My = masivtoch(Ny, qd, idy);
    else
        if HV % ���������� �������
            Mx = masivtoch(Nx, pd, idx);
            My = numHV;
        else % ���������� �����
            Mx = numHV;
            My = masivtoch(Ny, qd, idy);
        end
    end
else
    Mx = pixNx;
    My = pixNy;
end

% �������� ������� ����� �� ������� ����������
Nx = length(Mx);
Ny = length(My);

clear pd qd CR Npix pixNx pixNy HV numHV

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

clear vd wd Vd Wd f

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

% ����� ��� �� ������������� ����� ���� �����, ��������� ����� �������� ���� 
% �� ���� ������ (������ � �����)
            
%////////////////////
%// �������� ���� //
%////////////////////

% ���������� �������� �������� ���������� �������� ����
% ����������� ������� �� ������ �� ������ �� ������� ����
Rt = radzem(Rzmin, Rzmax, gama);

% ����������� ����� �������� ����
Rk = radkrivzem(Rzmin, Rzmax, gama);

clear gama

% ����������� ������� �� ������
dh = Rt - Rz;

clear Rt Rzmin Rzmax

% �������� ������
H = h + dh;

clear dh

% �������� ������ �� ����������� ������ �� ���� �������� �� ������� �������� ����
% ����� �� ������ �� ������� ���������� (������� ������� �� ����� ��������� ��)
Hrz = heighttoploschin(H, Rk, theta, phi);

% ���� �� ����� ���������, � ������ �� ������ � �������� (0�)
% �� ������� ������ (90�) � �������� ��������
deltaTheta = kytmizploschinami(0.0, theta, phi, H, Rk);
deltaPhi = kytmizploschinami(0.5 * pi, theta, phi, H, Rk);

% �������� ���� �������� ����� �� ����������� �������� ����
 theta_0 = theta + deltaTheta;
 phi_0 = phi + deltaPhi;
 
% ����������� ���������� ���������� ������ �������� ���
[l_theta, l_phi] = centerMPV(H, Rk, theta, phi, theta_0, phi_0);
 
% г����� ������ �� ������������ (�� ������� �������) � �������� �������� �� ���
deltaL_t = Hrz * tan(theta_0) - l_theta;
deltaL_p = Hrz * tan(phi_0) - l_phi;
 
% ���������� ���� ������� ��������� � ��������� �� ���� ������
for k = 1:5
    % �������� �� �������, � ��� �� �����
    [theta1{k}, phi1{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    [Theta1{k}, Phi1{k}] = naklonGlob(theta_0, phi_0, Theta0{k}, Phi0{k}, 0);
    
    % �������� �� �����, � ��� �� �������
    [theta2{k}, phi2{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    [Theta2{k}, Phi2{k}] = naklonGlob(theta_0, phi_0, Theta0{k}, Phi0{k}, 1);
    
    % �������� ������� ��� ���������� �������� naklonLok
    [theta3{k}, phi3{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    [Theta3{k}, Phi3{k}] = naklonLok(theta_0, phi_0, Theta0{k}, Phi0{k}, 0);
    
    % �������� ������� ��� ���������� ��������
    [theta4{k}, phi4{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    [Theta4{k}, Phi4{k}] = naklonLok(theta_0, phi_0, Theta0{k}, Phi0{k}, 1);
    
    if k < 5
        % �������� �� �������, � ��� �� �����
        [Rtheta1{k}, Rphi1{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
        % �������� �� �����, � ��� �� �������
        [Rtheta2{k}, Rphi2{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);
        
        % �������� ������� ��� ���������� ��������
        % �������� �� �������, � ��� �� �����
        [Rtheta3{k}, Rphi3{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
        % �������� �� �����, � ��� �� �������
        [Rtheta4{k}, Rphi4{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);
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
        
        % �������� ������� ��� ���������� ��������
        PhiR3(i) = acot((tan(theta3{1}(i, 1)) - tan(theta3{1}(i, Ny))) \...
            (tan(phi3{1}(i, 1)) - tan(phi3{1}(i, Ny))));
        % � ��������
        PhiR3(i) = 90 - PhiR3(i) * (180 / pi);
    else
        % � ��������
        PhiR1(i) = 90;
        PhiR3(i) = 90;
    end
    
    if (tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) ~= 0
        PhiR2(i) = acot((tan(theta2{1}(i, 1)) - tan(theta2{1}(i, Ny))) \...
            (tan(phi2{1}(i, 1)) - tan(phi2{1}(i, Ny))));
        % � ��������
        PhiR2(i) = 90 - PhiR2(i) * (180 / pi);
        
        % �������� ������� ��� ���������� ��������
        PhiR4(i) = acot((tan(theta4{1}(i, 1)) - tan(theta4{1}(i, Ny))) \...
            (tan(phi4{1}(i, 1)) - tan(phi4{1}(i, Ny))));
        % � ��������
        PhiR4(i) = 90 - PhiR4(i) * (180 / pi);
    else
        % � ��������
        PhiR2(i) = 90;
        PhiR4(i) = 90;
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
        
        % �������� ������� ��� ���������� ��������
        PhiC3(j) = atan((tan(theta3{1}(1, j)) - tan(theta3{1}(Nx, j))) \...
            (tan(phi3{1}(1, j)) - tan(phi3{1}(Nx, j))));
        % � ��������
        PhiC3(j) = PhiC3(j) * (180 / pi);
    else
        % � ��������
        PhiC1(j) = 0;
        PhiC3(j) = 0;
    end
    
    if (tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) ~= 0
        PhiC2(j) = atan((tan(theta2{1}(1, j)) - tan(theta2{1}(Nx, j))) \...
            (tan(phi2{1}(1, j)) - tan(phi2{1}(Nx, j))));
        % � ��������
        PhiC2(j) = PhiC2(j) * (180 / pi);
        
        % �������� ������� ��� ���������� ��������
        PhiC4(j) = atan((tan(theta4{1}(1, j)) - tan(theta4{1}(Nx, j))) \...
            (tan(phi4{1}(1, j)) - tan(phi4{1}(Nx, j))));
        % � ��������
        PhiC4(j) = PhiC4(j) * (180 / pi);
    else
        % � ��������
        PhiC2(j) = 0;
        PhiC4(j) = 0;
    end
end

clear j

% ����������� ���� ���������� �����
for k = 1:5
    % ��� �� ������
    Mtheta1{k} = Hrz .* tan(theta1{k}) + deltaL_t;
    Mtheta2{k} = Hrz .* tan(theta2{k}) + deltaL_t;
    
    Mphi1{k} = Hrz .* tan(phi1{k}) + deltaL_p;
    Mphi2{k} = Hrz .* tan(phi2{k}) + deltaL_p;
    
    % �������� ������� ��� ���������� ��������
    Mtheta3{k} = Hrz .* tan(theta3{k}) + deltaL_t;
    Mtheta4{k} = Hrz .* tan(theta4{k}) + deltaL_t;
    
    Mphi3{k} = Hrz .* tan(phi3{k}) + deltaL_p;
    Mphi4{k} = Hrz .* tan(phi4{k}) + deltaL_p;
    
    % ��� �� ������
    MTheta1{k} = Hrz .* tan(Theta1{k}) + deltaL_t;
    MTheta2{k} = Hrz .* tan(Theta2{k}) + deltaL_t;
    
    MPhi1{k} = Hrz .* tan(Phi1{k}) + deltaL_p;
    MPhi2{k} = Hrz .* tan(Phi2{k}) + deltaL_p;
    
    % �������� ������� ��� ���������� ��������
    MTheta3{k} = Hrz .* tan(Theta3{k}) + deltaL_t;
    MTheta4{k} = Hrz .* tan(Theta4{k}) + deltaL_t;
    
    MPhi3{k} = Hrz .* tan(Phi3{k}) + deltaL_p;
    MPhi4{k} = Hrz .* tan(Phi4{k}) + deltaL_p;
    
    % ��� �� �������� ��������
    if k < 5
        MRtheta1{k} = Hrz .* tan(Rtheta1{k}) + deltaL_t;
        MRtheta2{k} = Hrz .* tan(Rtheta2{k}) + deltaL_t;
        
        MRphi1{k} = Hrz .* tan(Rphi1{k}) + deltaL_p;
        MRphi2{k} = Hrz .* tan(Rphi2{k}) + deltaL_p;
        
        % �������� ������� ��� ���������� ��������
        MRtheta3{k} = Hrz .* tan(Rtheta3{k}) + deltaL_t;
        MRtheta4{k} = Hrz .* tan(Rtheta4{k}) + deltaL_t;
        
        MRphi3{k} = Hrz .* tan(Rphi3{k}) + deltaL_p;
        MRphi4{k} = Hrz .* tan(Rphi4{k}) + deltaL_p;
    end
end

clear k Theta1 Phi1 Theta2 Phi2 theta1 phi1 theta2 phi2...
    Theta3 Phi3 Theta4 Phi4 theta3 phi3 theta4 phi4

% ³���������� ����� ���������� ������, � �� ������ ������
if (id == 1) && (idx == 1) && (idy == 1) && (Nx0 > 1) && (Ny0 > 1)
    % �������� �ٲ�������
    % ���������� ��������� ����� � ����������� �����
    Mx = sysidtoch(Mx);
    My = sysidtoch(My);

    % ��������� ����� ��� ��� ������ (center weight - CW)
    [CWx1, CWy1] = centerweight(Mx, My, MTheta1, MPhi1, Nx, Ny);
    [CWx2, CWy2] = centerweight(Mx, My, MTheta2, MPhi2, Nx, Ny);
    
    % �������� ������� ��� ���������� ��������
    [CWx3, CWy3] = centerweight(Mx, My, MTheta3, MPhi3, Nx, Ny);
    [CWx4, CWy4] = centerweight(Mx, My, MTheta4, MPhi4, Nx, Ny);

    % ��������� ������� ��� ��������� ��������� ������� �����������
    [MinT1, MinP1] = minvidmizpix(Mx, My, CWx1, CWy1, Nx, Ny, proc);
    [MinT2, MinP2] = minvidmizpix(Mx, My, CWx2, CWy2, Nx, Ny, proc);
    
    % �������� ������� ��� ���������� ��������
    [MinT3, MinP3] = minvidmizpix(Mx, My, CWx3, CWy3, Nx, Ny, proc);
    [MinT4, MinP4] = minvidmizpix(Mx, My, CWx4, CWy4, Nx, Ny, proc);
    
    % ��������� ������������ ����. ��������, ��� ���� ��� ����� �������
    % �������������� �������� ������ � � ��������� ���� �� ������ ������
    
    % ������ - ����
    if ((isnan(MinT1) && ~isnan(MinP1)))
        MaxSc1 = MinP1;
    elseif (isnan(MinP1) && ~isnan(MinT1))
        MaxSc1 = MinT1;
    elseif ((isnan(MinT1) && isnan(MinP1)))
        MaxSc1 = 1.0;
    else
        MaxSc1 = max(MinT1, MinP1);
    end
    
    % ���� - ������
    if ((isnan(MinT2) && ~isnan(MinP2)))
        MaxSc2 = MinP2;
    elseif (isnan(MinP2) && ~isnan(MinT2))
        MaxSc2 = MinT2;
    elseif ((isnan(MinT2) && isnan(MinP2)))
        MaxSc2 = 1.0;
    else
        MaxSc2 = max(MinT2, MinP2);
    end
    
    % �������� ������� ��� ���������� ��������
    % ������ - ����
    if ((isnan(MinT3) && ~isnan(MinP3)))
        MaxSc3 = MinP3;
    elseif (isnan(MinP3) && ~isnan(MinT3))
        MaxSc3 = MinT3;
    elseif ((isnan(MinT3) && isnan(MinP3)))
        MaxSc3 = 1.0;
    else
        MaxSc3 = max(MinT3, MinP3);
    end
    
    % ���� - ������
    if ((isnan(MinT4) && ~isnan(MinP4)))
        MaxSc4 = MinP4;
    elseif (isnan(MinP4) && ~isnan(MinT4))
        MaxSc4 = MinT4;
    elseif ((isnan(MinT4) && isnan(MinP4)))
        MaxSc4 = 1.0;
    else
        MaxSc4 = max(MinT4, MinP4);
    end
    
    % ���������� �� ���������� ������������ ����� ������
    [Mtheta1, Mphi1] = everyscale(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    % �������� ������� ��� ���������� ��������
    [Mtheta3, Mphi3] = everyscale(Nx, Ny, MinT3, MinP3, Mtheta3, Mphi3,...
        CWx3, CWy3);
    [Mtheta4, Mphi4] = everyscale(Nx, Ny, MinT4, MinP4, Mtheta4, Mphi4,...
        CWx4, CWy4);
    
    % ��� ���������� ���� �� ������� 1 ��������������, ��� �������� ���
    % ����� ��� �������, ������� ��������� ��������� ��� ������ ���
    [CWx1, CWy1] = centerweight2(MTheta1, MPhi1, MinT1, MinP1,...
        Nx, Ny, CWx1, CWy1);
    [CWx2, CWy2] = centerweight2(MTheta2, MPhi2, MinT2, MinP2,...
        Nx, Ny, CWx2, CWy2);
    
    % �������� ������� ��� ���������� ��������
    [CWx3, CWy3] = centerweight2(MTheta3, MPhi3, MinT3, MinP3,...
        Nx, Ny, CWx3, CWy3);
    [CWx4, CWy4] = centerweight2(MTheta4, MPhi4, MinT4, MinP4,...
        Nx, Ny, CWx4, CWy4);
    
    % ��������� �������������� ������� ���, ��� ���������� �������������
    % �� ����� �� ����
    % ��������� ����������� ���� �������������
    Min1 = min(MinT1, MinP1);
    MinT1 = MinT1 / Min1;
    MinP1 = MinP1 / Min1;
    Min2 = min(MinT2, MinP2);
    MinT2 = MinT2 / Min2;
    MinP2 = MinP2 / Min2;
    
    % �������� ������� ��� ���������� ��������
    Min3 = min(MinT3, MinP3);
    MinT3 = MinT3 / Min3;
    MinP3 = MinP3 / Min3;
    Min4 = min(MinT4, MinP4);
    MinT4 = MinT4 / Min4;
    MinP4 = MinP4 / Min4;
    
    clear Min1 Min2 Min3 Min4
    
    % ���������� �� ���������� ������������ ����� ������
    [Mtheta1, Mphi1] = everyscale2(Nx, Ny, MinT1, MinP1, Mtheta1, Mphi1,...
        CWx1, CWy1);
    [Mtheta2, Mphi2] = everyscale2(Nx, Ny, MinT2, MinP2, Mtheta2, Mphi2,...
        CWx2, CWy2);
    
    % �������� ������� ��� ���������� ��������
    [Mtheta3, Mphi3] = everyscale2(Nx, Ny, MinT3, MinP3, Mtheta3, Mphi3,...
        CWx3, CWy3);
    [Mtheta4, Mphi4] = everyscale2(Nx, Ny, MinT4, MinP4, Mtheta4, Mphi4,...
        CWx4, CWy4);
    
    clear MinT1 MinP1 MinT2 MinP2 ...%CWx1 CWy1 CWx2 CWy2
        MinT3 MinP3 MinT4 MinP4
    
    % ���������� 3-� ��� � ������� ������� �����������
    [Mtheta1, Mphi1] = everyscale3(MaxSc1, Mtheta1, Mphi1);
    [Mtheta2, Mphi2] = everyscale3(MaxSc2, Mtheta2, Mphi2);
    
    % �������� ������� ��� ���������� ��������
    [Mtheta3, Mphi3] = everyscale3(MaxSc3, Mtheta3, Mphi3);
    [Mtheta4, Mphi4] = everyscale3(MaxSc4, Mtheta4, Mphi4);
end

clear idx idy Nx0 Ny0 Mx My

% ������� ���������� ��������
figure(1);
clf;
hold on;
grid off;   % �������� ���� "on", ��������� "off"
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (������ - ����) ����. + ���. �������');
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
    
    % �������� ������� ��� ���������� ��������
    ymin = min(min(Mphi3{k}(:)), ymin);
    ymax = max(max(Mphi3{k}(:)), ymax);
    xmin = min(min(Mtheta3{k}(:)), xmin);
    xmax = max(max(Mtheta3{k}(:)), xmax);
    
    ymin = min(min(Mphi4{k}(:)), ymin);
    ymax = max(max(Mphi4{k}(:)), ymax);
    xmin = min(min(Mtheta4{k}(:)), xmin);
    xmax = max(max(Mtheta4{k}(:)), xmax);
end

% ���������� ���������� ��
ylen = abs(ymax - ymin);
xlen = abs(xmax - xmin);

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

% figure(4);
% clf;
% hold on;
% grid on;   % �������� ���� "on", ��������� "off"
% xlabel('������� ������� [�]');
% ylabel('������ ������� [�]');
% title('�������� ������� ������ (���� - ������) ��������� �������');

% ylim([ymin ymax]);
% xlim([xmin xmax]);

% ylim(xcenter - 0.5 .* [len -len]);
% xlim(ycenter - 0.5 .* [len -len]);

% clear xcenter ycenter len

%  ³��������� �������� ������
for i = 1:Nx
   for j = 1:Ny
       plot(Mphi4{1}(i, j), Mtheta4{1}(i, j), 'bx');
       plot([Mphi4{2}(i, j), Mphi4{3}(i, j),...
           Mphi4{4}(i, j), Mphi4{5}(i, j), Mphi4{2}(i, j)],...
           [Mtheta4{2}(i, j), Mtheta4{3}(i, j), Mtheta4{4}(i, j),...
           Mtheta4{5}(i, j), Mtheta4{2}(i, j)], '-k');
       
       % �������� ���������� ���������� ����� ��� ��
       if id == 0
           for k = 1:4
               plot(MRphi4{k}(i, j), MRtheta4{k}(i, j), 'bx');
           end
       end
%        if id == 1
%            plot(CWy2{1}(i, j), CWx2{1}(i, j), 'b+');
%        end
   end
end

% �������� �� ������ ������ ������ �� ���� ��
% plot([Mphi4{1}(Nx, 1), Mphi4{1}(Nx, Ny),],...
%     [Mtheta4{1}(Nx, 1), Mtheta4{1}(Nx, Ny)], '-g');
% plot([Mphi4{1}(1, 1), Mphi4{1}(Nx, 1),],...
%     [Mtheta4{1}(1, 1), Mtheta4{1}(Nx, 1)], '-g');
% 
% plot([Mphi4{1}(1, 1), Mphi4{1}(1, Ny),],...
%     [Mtheta4{1}(1, 1), Mtheta4{1}(1, Ny)], '-g');

figure(2);
clf;
hold on;
grid off;   % �������� ���� "on", ��������� "off"
xlabel('������� ������� [�]');
ylabel('������ ������� [�]');
title('�������� ������� ������ (���� - ������) ����. + ���. �������');

% ylim([ymin ymax]);
% xlim([xmin xmax]);

ylim(xcenter - 0.5 .* [len -len]);
xlim(ycenter - 0.5 .* [len -len]);

% clear xcenter ycenter len

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

% �������� ������� ��� ���������� ��������
% figure(1);
% clf;
% hold on;
% grid on;   % �������� ���� "on", ��������� "off"
% xlabel('������� ������� [�]');
% ylabel('������ ������� [�]');
% title('�������� ������� ������ (������ - ����) ��������� �������');

% ylim([ymin ymax]);
% xlim([xmin xmax]);

% ylim(xcenter - 0.5 .* [len -len]);
% xlim(ycenter - 0.5 .* [len -len]);

% clear xcenter ycenter len

%  ³��������� �������� ������
for i = 1:Nx
   for j = 1:Ny
       plot(Mphi3{1}(i, j), Mtheta3{1}(i, j), 'bx');
       plot([Mphi3{2}(i, j), Mphi3{3}(i, j),...
           Mphi3{4}(i, j), Mphi3{5}(i, j), Mphi3{2}(i, j)],...
           [Mtheta3{2}(i, j), Mtheta3{3}(i, j), Mtheta3{4}(i, j),...
           Mtheta3{5}(i, j), Mtheta3{2}(i, j)], '-k');
       
       % �������� ���������� ���������� ����� ��� ��
       if id == 0
           for k = 1:4
               plot(MRphi3{k}(i, j), MRtheta3{k}(i, j), 'bx');
           end
       end
%        if id == 1
%            plot(CWy2{1}(i, j), CWx2{1}(i, j), 'b+');
%        end
   end
end

% �������� �� ������ ������ ������ �� ���� ��
% plot([Mphi3{1}(Nx, 1), Mphi3{1}(Nx, Ny),],...
%     [Mtheta3{1}(Nx, 1), Mtheta3{1}(Nx, Ny)], '-g');
% plot([Mphi3{1}(1, 1), Mphi3{1}(Nx, 1),],...
%     [Mtheta3{1}(1, 1), Mtheta3{1}(Nx, 1)], '-g');
% 
% plot([Mphi3{1}(1, 1), Mphi3{1}(1, Ny),],...
%     [Mtheta3{1}(1, 1), Mtheta3{1}(1, Ny)], '-g');

clear k i j pd qd id CWx1 CWy1 CWx2 CWy2 Mphi1 Mphi2 Mtheta1 Mtheta2...
    CWx3 CWy3 CWx4 CWy4 Mphi3 Mphi4 Mtheta3 Mtheta4

% ���������� ����� �� ������������� �������� �������� ������� (�� ������)
% 1-�� ������ (��� ���������� �������� � ������ �� ����)
% [Ltheta1, Lphi1] = razreshtoch(MTheta1, MPhi1, Nx, Ny);
% [Ltheta2, Lphi2] = razreshtoch(MTheta2, MPhi2, Nx, Ny);

% 2-�� ������ ���������� ������ �� �������� �������
% �������� �� �������, � ��� �� �����
Ltheta1 = vipstanmizpix(MRtheta1{1}, MRtheta1{3}, MRphi1{1}, MRphi1{3});
Lphi1 = vipstanmizpix(MRtheta1{2}, MRtheta1{4}, MRphi1{2}, MRphi1{4});

Ltheta3 = vipstanmizpix(MRtheta3{1}, MRtheta3{3}, MRphi3{1}, MRphi3{3});
Lphi3 = vipstanmizpix(MRtheta3{2}, MRtheta3{4}, MRphi3{2}, MRphi3{4});

% �������� �� �����, � ��� �� �������
Ltheta2 = vipstanmizpix(MRtheta2{1}, MRtheta2{3}, MRphi2{1}, MRphi2{3});
Lphi2 = vipstanmizpix(MRtheta2{2}, MRtheta2{4}, MRphi2{2}, MRphi2{4});

Ltheta4 = vipstanmizpix(MRtheta4{1}, MRtheta4{3}, MRphi4{1}, MRphi4{3});
Lphi4 = vipstanmizpix(MRtheta4{2}, MRtheta4{4}, MRphi4{2}, MRphi4{4});

clear MTheta3 MPhi3 MTheta4 MPhi4 MRtheta3 MRphi3 MRtheta4 MRphi4

% ���������� ������� [�� ^ -1]
NuTheta1 = 1000 ./ Ltheta1;
NuLphi1 = 1000 ./ Lphi1;
NuTheta2 = 1000 ./ Ltheta2;
NuLphi2 = 1000 ./ Lphi2;

NuTheta3 = 1000 ./ Ltheta3;
NuLphi3 = 1000 ./ Lphi3;
NuTheta4 = 1000 ./ Ltheta4;
NuLphi4 = 1000 ./ Lphi4;

% ��������� ������� 

clear NuTheta3 NuTheta4 NuLphi3 NuLphi4 PhiC3 PhiC4 PhiR3 PhiR4...
    NuTheta1 NuTheta2 NuLphi1 NuLphi2 PhiC1 PhiC2 PhiR1 PhiR2
    %Lphi1 Lphi2 Ltheta1 Ltheta2


% ����������� ���� ��������� �� ������
for i = 1:Nx
    for j = 1:Ny
        for k = 1:4
            alpha1{k}(i, j) = anglealpha(Rtheta1{k}(i, j), Rphi1{k}(i, j));
            alpha2{k}(i, j) = anglealpha(Rtheta2{k}(i, j), Rphi2{k}(i, j));
            
            alpha3{k}(i, j) = anglealpha(Rtheta3{k}(i, j), Rphi3{k}(i, j));
            alpha4{k}(i, j) = anglealpha(Rtheta4{k}(i, j), Rphi4{k}(i, j));
        end
    end
end

clear i j k Rtheta1 Rphi1 Rtheta2 Rphi2...
    Rtheta3 Rphi3 Rtheta4 Rphi4

% ��������� ������� �� �� �� ����� ����� ��� ����� �� ������
for k = 1:4
    dovzin1{k} = dovzdotoch(h, Rk, alpha1{k});
    dovzin2{k} = dovzdotoch(h, Rk, alpha2{k});
    
    dovzin3{k} = dovzdotoch(h, Rk, alpha3{k});
    dovzin4{k} = dovzdotoch(h, Rk, alpha4{k});
end

clear alpha1 alpha2 h k

% ��������� �������� �������� ���������� ��������, ��� ����������� ��
% �����
[LRatheta1, LRaphi1] = vidsmiztoch(dovzin1, Wx, Wy, Nx, Ny);
[LRatheta2, LRaphi2] = vidsmiztoch(dovzin2, Wx, Wy, Nx, Ny);

[LRatheta3, LRaphi3] = vidsmiztoch(dovzin3, Wx, Wy, Nx, Ny);
[LRatheta4, LRaphi4] = vidsmiztoch(dovzin4, Wx, Wy, Nx, Ny);

clear Nx Ny Wx Wy dovzin1 dovzin2 dovzin3 dovzin4

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

dlLRtheta1 = LRatheta1 - Ltheta1;
dlLRphi1 = LRaphi1 - Lphi1;
dlLRtheta2 = LRatheta2 - Ltheta2;
dlLRphi2 = LRaphi2 - Lphi2;

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

% roundn([min(dlLRtheta1(:)) max(dlLRtheta1(:)); min(dlLRphi1(:)) max(dlLRphi1(:))], -3)
% roundn([min(dlLRtheta2(:)) max(dlLRtheta2(:)); min(dlLRphi2(:)) max(dlLRphi2(:))], -3)

clear ans
