function [WX, WY] = angleRZ(Wx, Wy, Vd, Wd, f, Nx, Ny)
% ������� ���������� ������� ��������� ������

% WX, WY - ����� �������� ������� ������� ������ ��� ������ ������
% Wx, Wy - ����� ����� ������ ������
% Vd, Wd - ����� ������
% f - ������� �������

% �������� ���������� Wx_2, Wy_2
for i = 1:Nx
    for j = 1:Ny
        Wx_2(i, j) = atan(tan(Wy(j)) .* cos(Wx(i)));
        Wy_2(i, j) = atan(tan(Wx(i)) .* cos(Wy(j)));
    end
end

% ��� ����������� ����� �������� ������
for i = 1:Nx
    for j = 1:Ny
        WX(i, j) = acot((f ./ (Vd .* (cos(Wx_2(i, j)) .^ 2.0) .* cos(Wy_2(i, j))))...
            - ((Vd .* cos(Wy_2(i, j))) ./ (4.0 .* f)));
        WY(i, j) = acot((f ./ (Wd .* (cos(Wy_2(i, j)) .^ 2.0) .* cos(Wx_2(i, j))))...
            - ((Wd .* cos(Wx_2(i, j))) ./ (4.0 .* f)));
    end
end

end