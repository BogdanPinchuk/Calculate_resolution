function thetaZ = angleZemli(H, Rk, alpha)
% ������� ���������� ����������� ���� ����
% �� ���� �� ��������� � ����� ������ � ������������� � �����
% �������������

% H - ������� ������ ��������� � �����
% Rk - ����� �������� ���� � ��������� �����
% alpha ��� ��������� �� �� ������

% ��� �� ������ ������ � ������ �� �������� ����
thetaZ = asin(sin(alpha) * (H + Rk) / Rk) - alpha;

end