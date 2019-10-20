function [CWx, CWy] = centerweight(Kx, Ky, MTheta, MPhi, Nx, Ny)
% ������� ���������� ������ ��� ��� ������ "�����" �� ���

% Kx, Ky - ����� ������� ��������� ����� � ����������� �����
% MTheta, MPhi - ���� ���������� �����
% Nx, Ny - ������� ����� �� ������� ����������

% ��������� ������� �������� ���������� ����� �����

% ��������� min � max ������� � ����
for i = 1:Nx
    for j = 1:Ny
        % ���������� �� ������ ��� �����
        if ((Kx(1, i) - Kx(2, i)) == 0) && ((Ky(1, j) - Ky(2, j)) == 0)
            % �������� ������ ��� ���
%             CWx{1}(i, j) = MTheta{1}(i, j);
%             CWy{1}(i, j) = MPhi{1}(i, j);
            
            % �������� ����� �����
            for k = 1:5
                CWx{k}(i, j) = MTheta{k}(i, j);
                CWy{k}(i, j) = MPhi{k}(i, j);
            end
        else
            % �������� ����� �����
            CWx{2}(i, j) = MTheta{2}(Kx(2, i), Ky(1, j));
            CWy{2}(i, j) = MPhi{2}(Kx(2, i), Ky(1, j));
            
            CWx{3}(i, j) = MTheta{3}(Kx(2, i), Ky(2, j));
            CWy{3}(i, j) = MPhi{3}(Kx(2, i), Ky(2, j));
            
            CWx{4}(i, j) = MTheta{4}(Kx(1, i), Ky(2, j));
            CWy{4}(i, j) = MPhi{4}(Kx(1, i), Ky(2, j));
            
            CWx{5}(i, j) = MTheta{5}(Kx(1, i), Ky(1, j));
            CWy{5}(i, j) = MPhi{5}(Kx(1, i), Ky(1, j));
            
            % ���� � ������ � �����
            [CWx{1}(i, j), CWy{1}(i, j)] = peretintoch(CWx{2}(i, j),...
                CWy{2}(i, j), CWx{4}(i, j), CWy{4}(i, j), CWx{5}(i, j),...
                CWy{5}(i, j), CWx{3}(i, j), CWy{3}(i, j));
        end
    end
end

end
