function deltaAlpha = kytmizploschinami(fi, theta, phi, H, Rk)
% ������� ���������� ���� �� ����� ���������, 
% �� ����� ��� ������ ����� �� �������� �������

% fi - ��� �� ��������� ������� �� ��� �� ��� ��������� ������ ��� �� ����� ���������
% theta - ��� �������� � �������� �������
% phi - ��� �������� ������� �������

% ��������� ��� ��������� ��
alphaKA = anglealpha(theta, phi);

% ��������� ��������� ���� ��� � ����� ���� (��� �� ������ ������ � 
% ������ ��������� �� �� ������� ����)
thetaZ = angleZemli(H, Rk, alphaKA);

% ������ �� �������� ���������, ��������� ��� �� ��������� ��������������� ������� � �����������
% ��� (�� ����� ������ �� �� ���� �� ����� ������) �� ������ ��
if ((theta == 0) && ((0 <= phi) || (phi < 0.5 * pi)))
    phi_1 = 0;
elseif ((theta == 0) && ((-0.5 * pi < phi) || (phi < 0)))
    phi_1 = pi;
elseif ((theta ~= 0) && (abs(theta) < 0.5 * pi) && (phi == 0))
    phi_1 = sign(theta) * 0.5 * pi;
elseif (((0 < abs(theta)) || (abs(theta) < 0.5 * pi)) &&...
        ((0 < phi) || (phi < 0.5 * pi)))
    phi_1 = atan(tan(theta) / tan(phi));
elseif (((0 < abs(theta)) || (abs(theta) < 0.5 * pi)) &&...
        ((-0.5 * pi < phi) || (phi < 0)))
    phi_1 = pi + atan(tan(theta) / tan(phi));
else
    phi_1 = 0;
end

% ��� �� ����� ���������
deltaAlpha = 2 * asin(sin(phi_1 + fi) * sin(0.5 * thetaZ));

end