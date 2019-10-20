function [MinT, MinP] = minvidmizpix(Kx, Ky, CWx, CWy, Nx, Ny, proc)
% ������� ���������� ��������� ������� �� �������� ��������

% Kx, Ky - ����� ������� ��������� ����� � ����������� �����
% CWx, CWy - ������ ��� ��� "�����" ������ �� ���
% Nx, Ny - �� ������� ����������

% ������� ������� ���� �� ��������
% proc = 10;

% ������ �������� �������� ���������
% �� ��������
kx = 1;
ky = 1;
% �� ����������
tx = 1;
ty = 1;

Lx(1, 1) = NaN;
Ly(1, 1) = NaN;
MinT(1, 1) = NaN;
MinP(1, 1) = NaN;

% ��������� ������� �� ��������
for i = 1:Nx
    for j = 1:Ny
        % ���������� �� ������ ��� ����� ������ �� ��������
        % ((Kx(3, i) - Kx(3, i + 1)) ~= 1) - ����'������ ��������� ��������
        % � ����, ��� ������� �� ����� �� �������� ���� �������
        if (i < Nx) && (abs(Kx(3, i) - Kx(3, i + 1)) ~= 1)
            % ��������� ������� �� �������� ������� ������
            MinT(kx, tx) = vipstanmizpix(CWx{1}(i, j), CWx{1}(i + 1, j),...
                CWy{1}(i, j), CWy{1}(i + 1, j));
            
            % ��������� ������� ������ ������� �������� ������� ������
            % �� ������ ��� �'���� ������ ������
            [x, y] = peretintoch(CWx{2}(i, j), CWy{2}(i, j),...
                CWx{3}(i, j), CWy{3}(i, j), CWx{1}(i, j), CWy{1}(i, j),...
                CWx{1}(i + 1, j), CWy{1}(i + 1, j));
            
            % ��������� ������� �� ������� ����� �� ������ ������
            Lx(kx, tx) = vipstanmizpix(CWx{1}(i, j), x, CWy{1}(i, j), y);
            
            % ���������� ������� ������ ������� �������� ������� ������
            % �� ������ ��� �'���� ������ ������
            [x, y] = peretintoch(CWx{4}(i + 1, j), CWy{4}(i + 1, j),...
                CWx{5}(i + 1, j), CWy{5}(i + 1, j), CWx{1}(i, j),...
                CWy{1}(i, j), CWx{1}(i + 1, j), CWy{1}(i + 1, j));
            
            % ��������� ������� �� ������� ����� �� ������ ������
           Lx(kx, tx) = Lx(kx, tx) + vipstanmizpix(CWx{1}(i + 1, j),...
               x, CWy{1}(i + 1, j), y);
           
        end
        
        % ���������� �� ������ ��� ����� ������ �� ����������
        % ((Ky(3, j) - Ky(3, j + 1)) ~= 1) - ����'������ ��������� ��������
        % � ����, ��� ������� �� ����� �� �������� ���� �������
        if (j < Ny) && (abs(Ky(3, j) - Ky(3, j + 1)) ~= 1)
            % ��������� ������� �� �������� ��������
            MinP(ky, ty) = vipstanmizpix(CWx{1}(i, j), CWx{1}(i, j + 1),...
                CWy{1}(i, j), CWy{1}(i, j + 1));
            
            % ���������� ������� ������ ������� �������� ������� ������
            % �� ������ ��� �'���� ������ ������
            [x, y] = peretintoch(CWx{3}(i, j), CWy{3}(i, j),...
                CWx{4}(i, j), CWy{4}(i, j), CWx{1}(i, j), CWy{1}(i, j),...
                CWx{1}(i , j + 1), CWy{1}(i, j + 1));
            
            % ��������� ������� �� ������� ����� �� ������ ������
            Ly(ky, ty) = vipstanmizpix(CWx{1}(i, j), x, CWy{1}(i, j), y);
            
            % ���������� ������� ������ ������� �������� ������� ������
            % �� ������ ��� �'���� ������ ������
            [x, y] = peretintoch(CWx{2}(i, j + 1), CWy{2}(i, j + 1),...
                CWx{5}(i, j + 1), CWy{5}(i, j + 1), CWx{1}(i, j),...
                CWy{1}(i, j), CWx{1}(i, j + 1), CWy{1}(i, j + 1));
            
            % ��������� ������� �� ������� ����� �� ������ ������
            Ly(ky, ty) = Ly(ky, ty) + vipstanmizpix(CWx{1}(i, j + 1),...
                x, CWy{1}(i, j + 1), y);
        end
        
        % ��������� ��������� �� ����������
        if (j < Ny) && (abs((Ky(3, j) - Ky(3, j + 1))) ~= 1)
            tx = tx + 1;
            ty = ty + 1;
        else
            tx = tx + 1;
        end
    end
    
    % ��������� ��������� �� ��������
    if (i < Nx) && (abs((Kx(3, i) - Kx(3, i + 1))) ~= 1)
        kx = kx + 1;
        ky = ky + 1;
    else
        ky = ky + 1;
    end
    
    % ������� ��������
    tx = 1;
    ty = 1;
end

% ��������� ������� ������� ���� ���������� �� ��������
dLx = min(Lx(:)) * (proc / 100);
dLy = min(Ly(:)) * (proc / 100);
% dLx = Lx .* (proc / 100);
% dLy = Ly .* (proc / 100);

% ������ ������� �� ������� ��� �������� ����� �� �� ��������
Lx = Lx + dLx;
Ly = Ly + dLy;

% �������� ����������� �� �������� ������� �������������
MinT = MinT ./ Lx;
MinP = MinP ./ Ly;
% if isnan(Lx)
%     MinT = NaN;
% else
%     MinT = MinT ./ Lx;
% end
% if isnan(Ly)
%     MinP = NaN;
% else
%     MinP = MinP ./ Ly;
% end

% ��������� ��������� ����������
MinT = min(MinT(:));
MinP = min(MinP(:));

end
