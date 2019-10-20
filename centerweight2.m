function [CWx, CWy] = centerweight2(MTheta, MPhi, MinT, MinP,...
    Nx, Ny, CWx, CWy)
% ������� ���������� ������ ��� ��� ������ "�����" �� ���

% Mx, My - ����� ������� ��������� ����� � ����������� �����
% MTheta, MPhi - ���� ���������� �����
% MinT, MinP - �������� �����������
% Nx, Ny - ������� ����� �� ������� ����������

% % ������������� ������� �����
% [MTheta, MPhi] = everyscale(Nx, Ny, MinT, MinP, MTheta, MPhi,...
%     CWx, CWy);

% ��������� �������������� ������� ���, ��� ���������� �������������
% �� ����� �� ����
Min0 = min(MinT, MinP);
MinT = MinT / Min0;
MinP = MinP / Min0;

% ��������� �� ����� ��, �� �� 2-� �������������� ������ �������
% ���������� �������� ���� �� ����, ���� ��������������� ���������� ���
% ���������� ������ ���. ���� � ������ ������� ��� �������� ������
% ���������� � �������� ��������� ��������� �� ����� ����� ��� �������
% �������� ������ �� �������� ������. ���������� ����������� ����� ���
% �������� �������� ���� �� ������� ��� ������ ���

% �� ����� ���� �� � �����, �� ������ �� �������� ���, 
% ���� � ��������� �������������� ������ �����, � ��� ��� ��
% �� ���� ����������, ���� ��������� �������� ������ �� ����������
% � ��� ������ ��������� �� ����
cenX = round(Nx / 2);
cenY = round(Ny / 2);

% ���������� �� ��� �� ������ ���������� ������� �� �������� ������
% ������
if MinT == 1
    % ����������� ������ �� ������� ��� ������������ ��������
    S = vipstanmizpix(MTheta{1}(cenX, 1), MTheta{1}(cenX, Ny),...
        MPhi{1}(cenX, 1), MPhi{1}(cenX, Ny));
    k = cenX;
    
    % ����������� ������ �� ������� ��� ��� �������
    for i = 1:Nx
%         if i == 1
%             S = vipstanmizpix(MTheta{1}(i, 1), MTheta{1}(i, Ny),...
%                 MPhi{1}(i, 1), MPhi{1}(i, Ny));
%             k = i;
%         else
            if S > vipstanmizpix(MTheta{1}(i, 1), MTheta{1}(i, Ny),...
                    MPhi{1}(i, 1), MPhi{1}(i, Ny))
                S = vipstanmizpix(MTheta{1}(i, 1), MTheta{1}(i, Ny),...
                    MPhi{1}(i, 1), MPhi{1}(i, Ny));
                k = i;
            end
%         end
    end
else
    % ��������� ���� MinP == 1
    % ����������� ������ �� ������� ��� ������������ ��������
    S = vipstanmizpix(MTheta{1}(1, cenY), MTheta{1}(Nx, cenY),...
                MPhi{1}(1, cenY), MPhi{1}(Nx, cenY));
    k = cenY;
    
    % ����������� ������ �� ������� ��� ��� �������
    for j = 1:Ny
%         if j == 1
%             S = vipstanmizpix(MTheta{1}(1, j), MTheta{1}(Nx, j),...
%                 MPhi{1}(1, j), MPhi{1}(Nx, j));
%             k = j;
%         else
            if S > vipstanmizpix(MTheta{1}(1, j), MTheta{1}(Nx, j),...
                    MPhi{1}(1, j), MPhi{1}(Nx, j))
                S = vipstanmizpix(MTheta{1}(1, j), MTheta{1}(Nx, j),...
                    MPhi{1}(1, j), MPhi{1}(Nx, j));
                k = j;
            end
%         end
    end
end

% ��������� min � max ������� � ����
for i = 1:Nx
    for j = 1:Ny
        % ���������� �� ��� �� ������ "��������"
        if MinT == 1
            % �������� ������ ���
            CWx{1}(i, j) = CWx{1}(k, j);
            CWy{1}(i, j) = CWy{1}(k, j);
        else
            % ��������� ���� MinP == 1
            % �������� ������ ���
            CWx{1}(i, j) = CWx{1}(i, k);
            CWy{1}(i, j) = CWy{1}(i, k);
        end
    end
end

end
