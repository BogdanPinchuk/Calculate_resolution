function [LRatheta, LRaphi] = vidsmiztoch(dovzin, Wx, Wy, Nx, Ny)
% ������� ���������� �������� �� ������� ���������� ��������,
% ��� ������������� �����

% LRatheta - �� ������ ������� (������ �������� �������)
% LRaphi - �� ������ ����� (������� �������� �������)
% dovzin - ������� �� �� �� ����� ����� ��� ����� �� ������
% Wx, Wy - ����� �������� ������ (�������)

% ���� ��������� � ����� ���������� �� �������� �������� ��
for i = 1:Nx
    for j = 1:Ny
        LRatheta(i, j) = sqrt(dovzin{1}(i, j) .^ 2 + ...
            dovzin{3}(i, j) .^ 2 - 2 .* dovzin{1}(i, j) .*...
            dovzin{3}(i, j) .* cos(Wx(i)));
        
        LRaphi(i, j) = sqrt(dovzin{2}(i, j) .^ 2 + ...
            dovzin{4}(i, j) .^ 2 - 2 .* dovzin{2}(i, j) .*...
            dovzin{4}(i, j) .* cos(Wy(j)));
    end

end
