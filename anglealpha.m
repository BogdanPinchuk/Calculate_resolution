function alpha = anglealpha(theta, phi)
% ������� ���� ��������� �� ������

% alpha - ��� ��������� ����� �� �� ������
% theta - �������� ���� ��������� �� ������� ������ �������
% phi - �������� ���� ��������� �� ������� ������� �������
% �� ���� � �������

% ������ ������� ����� ������������ ��� ��� � 90�
% alpha = atan(sqrt(tan(theta) .^ 2 + tan(phi) .^ 2));

% �������� ������� ����������
% theta2 = atan(tan(phi) .* cos(theta));
% phi2 = atan(tan(theta) .* cos(phi));

% alpha = acos(cos(theta) .* cos(theta2));
% alpha = acos(cos(phi) .* cos(phi2));
% alpha = asin(sqrt(sin(theta2) .^ 2 + sin(phi2) .^ 2));

% ���� �������� ��������
if abs(phi) == (pi / 2)
    % ���� ��� phi = 90�, ��
    phi2 = atan(tan(theta) .* cos(phi));
    alpha = acos(cos(phi) .* cos(phi2));
elseif abs(theta) == (pi / 2)
    % ���� ��� theta = 90�, ��
    theta2 = atan(tan(phi) .* cos(theta));
    alpha = acos(cos(theta) .* cos(theta2));
elseif (abs(phi) ==(pi / 2)) && (abs(theta) == (pi / 2))
    alpha = pi / 2;
else
    alpha = atan(sqrt(tan(theta) .^ 2 + tan(phi) .^ 2));
end

end
