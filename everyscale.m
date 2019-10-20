function [MTheta, MPhi] = everyscale(Nx, Ny, MinT, MinP, Mtheta, Mphi,...
    CWx, CWy)
% ������� �������� ������������� ������

% Mtheta, Mphi - ���������� ����� �������� ������
% MinT, MinP - �������� �����������
% Nx, Ny - ������� ����� �� ������� ����������
% CWx, CWy - ������ ��� ��� "�����" ������ �� ���

% ��������� ������ �� ��������
if ~isnan(MinT) && ~isnan(MinP)
    MinSc = min(MinT, MinP);
else
    MinSc = 1.0;
end

% ����������
if MinSc ~= 1.0
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
else
    MTheta = Mtheta;
    MPhi = Mphi;
end
    
end
