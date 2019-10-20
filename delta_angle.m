function masiv = delta_angle(angle0, angle1, Nt)
% ������� ���������� ������ ��

% angle0 - ���������� ��� ����� ����������
% angle1 - �������� ��� ���� ���������
% Nt - ������� ����� �� �� ��������� ������� ���������
% masiv - ����� ������������ �����

% �������� �� ������� �����
if Nt < 2
    Nt = 2;
end

% ��������� ����
delta = (angle1 - angle0) ./ (Nt - 1);

for i = 1:Nt
    masiv(i) = angle0 + delta .* (i - 1);
end

end
