script % Projection_of_pixels_7
% ������ ������ F5
% ������ ��������� F9

%%%%%%%%
% ���� %
%%%%%%%%
clear; clc;

% ʳ������ �������� �������� ��������
% � �������� ������� (�� ��������, Ox) TDI
pd = 128;

% ������� ������� (�� ����������, Oy)
qd = 24001;

% ����� ������

% ����� ��������� ��������
% ������� [���]
vd = 7 * 10 ^ -6;  %8

% ������ [���]
wd = 7 * 10 ^ -6;  %8

% ����� �������� ��������
% ������ [���]
Vd = 7 * 10 ^ -6;  %8
% �������� [���]
Wd = 7 * 10 ^ -6;  %8

% ������� ������� ��'������ [��]
f = 6680  * 10 ^ -3; %2260

% ������� ������ ���� ��������� � ����� �� (��� � �����)
% ������ [����]
dAp = 0 * (pi / 180.0);
% �������� [����]
dAq = 0 * (pi / 180.0);

% ������ ����� [��]
h = 668 * 10 ^ 3;

% ���� ��������
% ������� [����]
Theta = 35 * (pi / 180.0);
% ����� [����]
Phi = 35 * (pi / 180.0);
% �������� [����]
psi = 0 * (pi / 180.0);

% ��������� ��� ������ ��������� ������� ������
% ʳ������ �����
Ntk = 8;

% ������ (-90�...+90�) �� ��� ����������� �� [����]
gama = 45 * (pi / 180.0);

% ������ ����, ����� ������� �� ����� ������� ���� �� �� ������ [��]
Rz = 6371.032 * 10 ^ 3;
Rzmin = 6356.777 * 10 ^ 3;
Rzmax = 6378.160 * 10 ^ 3;

%%%%%%%%%%%%%%%%%%%%%%%%
% ��������� ��������� %
%%%%%%%%%%%%%%%%%%%%%%%%

% ʳ������ ������� ������ ��� �����������
% ������ �������� �������
Nx0 = 3;

% ������� �������� �������
Ny0 = 3;

% ³���������� �����
% ���� id = 0 - ���������� ��
% ���� id = 1 - ���������� ����� ������
idx = 1;
idy = 1;
% ���������� ������������� (0 - �, 1 - ���)
id = 1;

% ������� �� �������� � % (��������)
proc = 10;

% ³���������� �������
graph = false;

% ������� ���������� 1 (false) ��� 
% ����� ������ (true)
Npix = true;

% ���������� ��� ��� true ������� �� ����� false
CR = true;

% ���������� ������� true, ��� ����� false
HV = false;

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

% clear dAp dAq

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

clear CR Npix pixNx pixNy HV numHV

% ���������� ������ ��������� ������ ������
for i = 1:Nx
%     Lpi(i) = Lp0 + (Mx(i) - 1) * Vd;
    
    % ����������� �����, ��� ����������� ��� ������� ������������ �� ��
    Lpi(i) = Lp0 + (Mx(Nx + 1 - i) - 1) * Vd;
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

% clear vd wd

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

clear mx my wx wy dwx dwy dWx dWy Wx Wy

% � ��������� �������� ����� ����� ��� ����� ������
% ������� ����� ���� ��������
for k = 1:5
    for i = 1:Nx
       for j = 1:Ny
           % ��� �� ������
%            theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
%                - tan(Wyi{k}(1, j)) * sin(psi));
%            phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
%                + tan(Wxi{k}(1, i)) * sin(psi));
%            
%            % ��� �񳺿 ������� (��) ������
%            Theta0{k}(i, j) = atan(tan(WXi{k}(1, i)) * cos(psi)...
%                - tan(WYi{k}(1, j)) * sin(psi));
%            Phi0{k}(i, j) = atan(tan(WYi{k}(1, j)) * cos(psi)...
%                + tan(WXi{k}(1, i)) * sin(psi));
           
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

clear k i j Wxi Wyi wxi wyi %WXi WYi psi

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

% clear gama

% ����������� ������� �� ������
dh = Rt - Rz;

clear Rt Rzmin Rzmax Rz

% �������� ������
H = h + dh;

clear dh %h

%%%%%%%%%%%%%%%%%%%%%%%%%
% �������
% ����� ���������� ������ ������� �������� �� �������� � ������
%%%%%%%%%%%%%%%%%%%%%%%%%

% ���� �������
thetaE = Theta / (Ntk - 1);
% ���� �����
phiE = Phi / (Ntk - 1);

for it = 1:Ntk
    for ik = 1:Ntk
        theta = (it - 1) * thetaE;
        phi = (ik - 1) * phiE;
        
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

        clear deltaPhi deltaTheta

        % ����������� ���������� ���������� ������ �������� ���
        [l_theta, l_phi] = centerMPV(H, Rk, theta, phi, theta_0, phi_0);

%         clear H Rk

        % г����� ������ �� ������������ (�� ������� �������) � �������� �������� �� ���
        deltaL_t = Hrz * tan(theta_0) - l_theta;
        deltaL_p = Hrz * tan(phi_0) - l_phi;

        clear l_phi l_theta

        % ���������� ���� ������� ��������� � ��������� �� ���� ������
        for k = 1:5
            if k < 5
                % �������� �� �������, � ��� �� �����
        %         [Rtheta1{k}, Rphi1{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
                % �������� �� �����, � ��� �� �������
        %         [Rtheta2{k}, Rphi2{k}] = naklonGlob(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);

                % �������� ������� ��� ���������� ��������
                % �������� �� �������, � ��� �� �����
                [Rtheta3{k}, Rphi3{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 0);
                % �������� �� �����, � ��� �� �������
%                 [Rtheta4{k}, Rphi4{k}] = naklonLok(theta_0, phi_0, Rtheta0{k}, Rphi0{k}, 1);
            end
        end

        clear k %theta0 phi0 Rtheta0 Rphi0
            %phi_0 theta_0 Theta0 Phi0 theta phi

        % ����������� ���� ���������� �����
        for k = 1:5
            % ��� �� �������� ��������
            if k < 5
        %         MRtheta1{k} = Hrz .* tan(Rtheta1{k}) + deltaL_t;
        %         MRtheta2{k} = Hrz .* tan(Rtheta2{k}) + deltaL_t;
        %         
        %         MRphi1{k} = Hrz .* tan(Rphi1{k}) + deltaL_p;
        %         MRphi2{k} = Hrz .* tan(Rphi2{k}) + deltaL_p;

                % �������� ������� ��� ���������� ��������
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

        % ���������� ����� �� ������������� �������� �������� ������� (�� ������)
        % �������� �� �������, � ��� �� �����
        % Ltheta1 = vipstanmizpix(MRtheta1{1}, MRtheta1{3}, MRphi1{1}, MRphi1{3});
        % Lphi1 = vipstanmizpix(MRtheta1{2}, MRtheta1{4}, MRphi1{2}, MRphi1{4});

        Ltheta3 = vipstanmizpix(MRtheta3{1}, MRtheta3{3}, MRphi3{1}, MRphi3{3});
        Lphi3 = vipstanmizpix(MRtheta3{2}, MRtheta3{4}, MRphi3{2}, MRphi3{4});

        % �������� �� �����, � ��� �� �������
        % Ltheta2 = vipstanmizpix(MRtheta2{1}, MRtheta2{3}, MRphi2{1}, MRphi2{3});
        % Lphi2 = vipstanmizpix(MRtheta2{2}, MRtheta2{4}, MRphi2{2}, MRphi2{4});

%         Ltheta4 = vipstanmizpix(MRtheta4{1}, MRtheta4{3}, MRphi4{1}, MRphi4{3});
%         Lphi4 = vipstanmizpix(MRtheta4{2}, MRtheta4{4}, MRphi4{2}, MRphi4{4});
        
        clear MTheta3 MPhi3 MTheta4 MPhi4 MRtheta3 MRphi3 MRtheta4 MRphi4...
            MRphi1 MRphi2 MRtheta1 MRtheta2

        % ���������� ������������ ��������� �� ���������� ���������� ����� �.�.
        % � �������� �������
        A = Vd .* Hrz ./ f;
        % ������� �������
        P = Wd .* Hrz ./ f;

        clear Hrz %Vd Wd f

        % �������� ���� ��������� ��� ��������
        theta_2 = atan(tan(phi_0) .* cos(theta_0));
        phi_2 = atan(tan(theta_0) .* cos(phi_0));

        for i = 1:Nx
           for j = 1:Ny
               % ��� �������� ��������

               %%%%%%%%%%%%%%%%%%%%%
               % ������� ������-���� (���. ���.)
               % ������� ����-������ (����. ���.)
               %%%%%%%%%%%%%%%%%%%%%

               % � �������� �������
        %        A3(i, j) = A ./ (cos(theta_0 + Theta0{1}(i, j)) .^ 2 .* ...
        %            cos(theta_2 + Phi0{1}(i, j)));
               A3(i, j) = A ./ (cos(theta_0 + WXi{1}(i)) .^ 2 .* ...
                   cos(theta_2 + WYi{1}(j)));
               % ������� �������
        %        P3(i, j) = P ./ (cos(theta_0 + Theta0{1}(i, j)) .* ...
        %            cos(theta_2 + Phi0{1}(i, j)) .^ 2);
               P3(i, j) = P ./ (cos(theta_0 + WXi{1}(i)) .* ...
                   cos(theta_2 + WYi{1}(j)) .^ 2);

               %%%%%%%%%%%%%%%%%%%%%
               % ������� ����-������ (���. ���.)
               % ������� ������-���� (����. ���.)
               %%%%%%%%%%%%%%%%%%%%%

               % � �������� �������
        %        A4(i, j) = A ./ (cos(phi_2 + Theta0{1}(i, j)) .^2 .* ...
        %            cos(phi_0 + Phi0{1}(i, j)));
%                A4(i, j) = A ./ (cos(phi_2 + WXi{1}(i)) .^2 .* ...
%                    cos(phi_0 + WYi{1}(j)));
               % ������� �������
        %        P4(i, j) = P ./ (cos(phi_2 + Theta0{1}(i, j)) .* ...
        %            cos(phi_0 + Phi0{1}(i, j)) .^2);
%                P4(i, j) = P ./ (cos(phi_2 + WXi{1}(i)) .* ...
%                    cos(phi_0 + WYi{1}(j)) .^2);
           end
        end

        clear theta_2 phi_2 theta_0 phi_0 A P Theta0 Phi0 %WXi WYi
        clear idx idy Nx0 Ny0 Mx My
        
        % ��������� ���������        
        dAbs3t = abs(Ltheta3 - A3);
        dAbs3k = abs(Lphi3 - P3);
%         dAbs4t = abs(Ltheta4 - A4);
%         dAbs4k = abs(Lphi4 - P4);
        
        % ³������ ���������
        dVid3t = 100 .* abs(1.0 - (A3 ./ Ltheta3));
        dVid3k = 100 .* abs(1.0 - (P3 ./ Lphi3));
%         dVid4t = 100 .* (1.0 - abs(A4 ./ Ltheta4));
%         dVid4k = 100 .* (1.0 - abs(P4 ./ Lphi4));
        
        % �������� � ��������� ����� ��� ���������
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

% ���������������� ����� ��� �����������
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
% ʳ����
% ����� ���������� ������ ������� �������� �� �������� � ������
%%%%%%%%%%%%%%%%%%%%%%%%%

clear k i j id CWx1 CWy1 CWx2 CWy2 Mphi1 Mphi2 Mtheta1 Mtheta2...
    CWx3 CWy3 CWx4 CWy4 Mphi3 Mphi4 Mtheta3 Mtheta4...
    len xcenter ycenter
clc;
clear ans graph proc

% �����������
% ����� �����
filename = [strrep(['��� �� ���������� �������� ', datestr(now)],...
    ':', '.'), '.xlsx'];

% filename = 'Test.xlsx';

% ������� ��� ����� ������� ��� ������� ���������
% ����� ��� ����
zvit = {};
zvit{size(zvit, 1) + 1, 2} = strrep(['��� �� ���������� �������� ',...
    datestr(now)], ':', '.');
zvit{size(zvit, 1) + 1, 2} = '������� ��������� ��� ��� - ������-����';
zvit{size(zvit, 1) + 1, 1} = [];

zvit{size(zvit, 1) + 1, 1} = pd;
zvit{size(zvit, 1), 2} = 'ʳ������ �������� �������� �������� �� ��������';
zvit{size(zvit, 1) + 1, 1} = qd;
zvit{size(zvit, 1), 2} = 'ʳ������ �������� �������� �������� �� ����������';
zvit{size(zvit, 1) + 1, 1} = vd * 10^6;
zvit{size(zvit, 1), 2} = '����� ��������� ��������, ������� [���]';
zvit{size(zvit, 1) + 1, 1} = wd * 10^6;
zvit{size(zvit, 1), 2} = '����� ��������� ��������, ������ [���]';
zvit{size(zvit, 1) + 1, 1} = Vd * 10^6;
zvit{size(zvit, 1), 2} = '����� �������� ��������, ������ [���]';
zvit{size(zvit, 1) + 1, 1} = Wd * 10^6;
zvit{size(zvit, 1), 2} = '����� ��������� ��������, �������� [���]';
zvit{size(zvit, 1) + 1, 1} = f * 10^3;
zvit{size(zvit, 1), 2} = '������� ������� ��`������ [��]';
zvit{size(zvit, 1) + 1, 1} = dAp * (180.0 / pi);
zvit{size(zvit, 1), 2} = '������� ������ ���� ���������, ������ [����]';
zvit{size(zvit, 1) + 1, 1} = dAq * (180.0 / pi);
zvit{size(zvit, 1), 2} = '������� ������ ���� ���������, �������� [����]';
zvit{size(zvit, 1) + 1, 1} = h * 10^(-3);
zvit{size(zvit, 1), 2} = '������ ����� [��]';
zvit{size(zvit, 1) + 1, 1} = Theta * (180.0 / pi);
zvit{size(zvit, 1), 2} = '��� ��������, ������ [����]';
zvit{size(zvit, 1) + 1, 1} = Phi * (180.0 / pi);
zvit{size(zvit, 1), 2} = '��� ��������, ���� [����]';
zvit{size(zvit, 1) + 1, 1} = psi * (180.0 / pi);
zvit{size(zvit, 1), 2} = '��� ��������, �������� [����]';
zvit{size(zvit, 1) + 1, 1} = gama * (180.0 / pi);
zvit{size(zvit, 1), 2} = '������ ���� �� ��� ����������� �� [����]';
zvit{size(zvit, 1) + 1, 1} = [];

zvit{size(zvit, 1) + 1, 2} = '������� ������ �������� �� �������� �� ����';
zvit{size(zvit, 1) + 1, 2} = '� �������� ������� 1..Np ����� �� �����';
zvit{size(zvit, 1) + 1, 2} = '������� ������� 1..Nq ���� �� �����';

clear pd qd vd wd Vd Wd f dAp dAq h theta phi psi gama

% �������� � ����
xlswrite(filename, zvit);

% ���������� �������� �� ����� � �������� ����� ������ �� ��������
% ������
for i = 1:Nx
   for j = 1:Ny
       % ������� ��� ��������� �������� �������
       zvit = {};
       
       y = 1;
       x = 1;
       
       % ³����� �������� � �������� �������
       zvit{x, y} =...
           '³������ ��������� � �������� �������, [%]';

       x = 2;
       
       % �������� ���
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dvid3t{i, j}(it, ik);
          end
       end
       
       y = Ntk + 1 + 2;
       x = 1;
       
       % ³����� �������� � ������� �������
       zvit{x, y} =...
           '³������ ��������� � ������� �������, [%]';
       
       x = 2;
       
       % �������� ���
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dvid3k{i, j}(it, ik);
          end
       end
       
       y = 1;
       x = Ntk + 1 + 3;
       
       % �������� �������� � �������� �������
       zvit{x, y} =...
           '��������� ��������� � �������� �������, [�]';
       
       x = Ntk + 1 + 4;
       
       % �������� ���
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dabs3t{i, j}(it, ik);
          end
       end
       
       y = Ntk + 1 + 2;
       x = Ntk + 1 + 3;
       
       % ³����� �������� � ������� �������
       zvit{x, y} =...
           '��������� ��������� � ������� �������, [�]';
       
       x = Ntk + 1 + 4;
       
       % �������� ���
       for it = 1:Ntk + 1
          for ik = 1:Ntk + 1
              zvit{it + x - 1, ik + y - 1} = Dabs3k{i, j}(it, ik);
          end
       end
       
       % ����� � ����
       sheetname = strcat('ϳ����� (', string(i), ', ', string(j), ')');
       xlswrite(filename, zvit, sheetname);
   end
end

clear Dabs3k Dabs3t Dvid3k Dvid3t

% �������� ���
% ������� ��� ��������� �������� �������
zvit = {};

for it = 1:Ntk
   for ik = 1:Ntk
       theta = round((it - 1) * thetaE * (180.0 / pi), 5);
       phi = round((ik - 1) * phiE * (180.0 / pi), 5);
       
       y = 1;
       x = 1;
       
       % ³����� �������� � �������� �������
       zvit{x, y} =...
           '������ ����. ������� ��������� � �������� �������, [%]';

       x = 2;
       
       % �������� ���
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaVid3t{it, ik}(:)) / (Nx * Ny);
       
       y = Ntk + 1 + 2;
       x = 1;
       
       % ³����� �������� � ������� �������
       zvit{x, y} =...
           '������ ����. ������� ��������� � ������� �������, [%]';
       
       x = 2;
       
       % �������� ���
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaVid3k{it, ik}(:)) / (Nx * Ny);
       
       y = 1;
       x = Ntk + 1 + 3;
       
       % �������� �������� � �������� �������
       zvit{x, y} =...
           '������ ����. ��������� ��������� � �������� �������, [�]';
       
       x = Ntk + 1 + 4;
       
       % �������� ���
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaAbs3t{it, ik}(:)) / (Nx * Ny);
       
       y = Ntk + 1 + 2;
       x = Ntk + 1 + 3;
       
       % ³����� �������� � ������� �������
       zvit{x, y} =...
           '������ ����. ��������� ��������� � ������� �������, [�]';
       
       x = Ntk + 1 + 4;
       
       % �������� ���
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sum(deltaAbs3k{it, ik}(:)) / (Nx * Ny);
       
       % �������� ������� ���������
       y = 2 * (Ntk + 1 + 2) - 1;
       x = 1;
       
       % ³����� �������� � �������� �������
       zvit{x, y} =...
           '������ ����. ������� ���������, [%]';
       
       x = 2;
       
       % �������� ���
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
       
       % ³����� �������� � ������� �������
       zvit{x, y} =...
           '������ ����. ��������� ���������, [�]';
       
       x = Ntk + 1 + 4;
       
       % �������� ���
       zvit{it + x, y} = theta;
       zvit{x, ik + y} = phi;
       zvit{it + x, ik + y} = ...
           sqrt(sum(deltaAbs3t{it, ik}(:)) *...
           sum(deltaAbs3k{it, ik}(:))) / (Nx * Ny);
       
   end
end

clear deltaAbs3k deltaAbs3t deltaVid3k deltaVid3t phi theta phiE thetaE

% ����� � ����
xlswrite(filename, zvit, '�����');

clear Ntk i j x y ik it

% �������� ��������� ��� ���������� ����
msgbox({'���������� �������� ���������� � ����:'...
    filename}, '���', 'help');

clear filename Phi Theta Nx Ny sheetname zvit
