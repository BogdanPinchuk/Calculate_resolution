function dovzin = dovzdotoch(h, Rk, alpha)
% ������� ���������� ������ �� �� �� ����� ������ �� ����

% dovzin - ������� �� �� �� ����� ����� ��� ����� �� ������
% h - ����������� ������ ��������� 
% Rk - ����� �������� ����
% alpha - ��� ��������� �� ������

beta = asin(((h + Rk) ./ Rk) .* sin(alpha)) - alpha;

dovzin = (h + Rk .* (1 - cos(beta))) ./ cos(alpha);

end
