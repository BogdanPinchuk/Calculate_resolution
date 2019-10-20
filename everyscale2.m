function [MTheta, MPhi] = everyscale2(Nx, Ny, MinT, MinP, Mtheta, Mphi,...
    CWx, CWy)
% ������� ���������� ������ ��� ��� ������ "�����" �� ���

% MTheta, MPhi - ��������������� ���������� ����� ������
% Mtheta, Mphi - ���������� ����� ������
% MinT, MinP - �������� �����������
% Nx, Ny - ������� ����� �� ������� ����������
% CWx, CWy - ������ ��� ��� "�����" ������ �� ���

% ��������� ���������� �������
if MinT == 1
    % ��������� ������� �� �� OY
    % ���� MinP = NaN (�������������)
    if isnan(MinP)
        MinSc = 1.0;
    else
        MinSc = MinP;
    end
else
    % ��������� ������� �� �� OX (�� MinP == 1)
    % ���� MinT = NaN (�������������)
    if isnan(MinT)
        MinSc = 1.0;
    else
        MinSc = MinT;
    end
end

% ����������
for i = 1:Nx
    for j = 1:Ny
        for k = 1:5
            MTheta{k}(i, j) = (MinSc * (Mtheta{k}(i, j) - CWx{1}(i, j)))...
                + CWx{1}(i, j);
            MPhi{k}(i, j) = (MinSc * (Mphi{k}(i, j) - CWy{1}(i, j)))...
                + CWy{1}(i, j);
        end
    end
end
    
end
