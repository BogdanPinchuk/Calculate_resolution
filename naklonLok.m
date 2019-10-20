function [Theta, Phi] = naklonLok(theta, phi, Wx, Wy, ind)
% ������� ���������� ������ ��

% Theta - theta1
% Phi - phi1
% theta, phi - ��������
% Wx - theta0
% Wy - phi0
% ind - ������ ���� ������� �� ������������ �������
% ���� 0 - �� �������� �� �������, � ��� �� �����
% ���� 1 - �� �������� �� �����, � ��� �� �������

% ��������� ����������
theta_2 = atan(tan(phi) .* cos(theta));
phi_2 = atan(tan(theta) .* cos(phi));

Wx_2 = atan(tan(Wy) .* cos(Wx));
Wy_2 = atan(tan(Wx) .* cos(Wy));

if ind == 0
   % ���� ind = 0
   % ��� ���������� ��������
   dtheta = atan(tan(Wy_2) ./ cos(theta_2 + Wy));
   Theta = theta + dtheta;
   Phi = atan(tan(theta_2 + Wy) .* cos(dtheta) ./ cos(Theta));
else
   % ���� ind = 1
   % ��� ���������� ��������
   dphi = atan(tan(Wx_2) ./ cos(phi_2 + Wx));
   Phi = phi + dphi;
   Theta = atan(tan(phi_2 + Wx) .* cos(dphi) ./ cos(Phi));
end

Theta = real(Theta);
Phi = real(Phi);

end
