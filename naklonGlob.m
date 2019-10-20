function [Theta, Phi] = naklonGlob(theta, phi, Wx, Wy, ind)
% ������� ���������� ������ ��

% Theta - theta1
% Phi - phi1
% theta, phi - ��������
% Wx - theta0
% Wy - phi0
% id - ������ ���� ������� �� ������������ �������
% ���� 0 - �� �������� �� �������, � ��� �� �����
% ���� 1 - �� �������� �� �����, � ��� �� �������

% ��������� ����������
theta_2 = atan(tan(phi) .* cos(theta));
phi_2 = atan(tan(theta) .* cos(phi));

if ind == 0
   % ���� ind = 0
   dphi = atan((tan(Wy) .* cos(Wx)) ./ cos(Wx + phi_2));
   Phi = phi + dphi;
   Theta = atan((tan(Wx + phi_2) .* cos(dphi)) ./ cos(Phi));
else
   % ���� ind = 1
   dtheta = atan((tan(Wx) .* cos(Wy)) ./ cos(Wy + theta_2));
   Theta = theta + dtheta;
   Phi = atan((tan(Wy + theta_2) .* cos(dtheta)) ./ cos(Theta));
end

Theta = real(Theta);
Phi = real(Phi);

end
