function Cz = masivtoch(nz, Nz, id)
% ������� ���������� ������� ��������� ������

% Nz - ������� ������ ������� �� ��������� ���������
% nz - ������� ������ ���������� ������������
% Cz - ����� �������� ����� ��� �����������
% id - ������������� ����������� �����
% ���� id = 0 - �� �������������
% ���� id = 1 - �������������

if id == 1
    % ����������� ���������
    for i = 1:nz
        if i == 1
            % ���� ��� ������� ���� ����� (����������)
            if (nz == 1) && (Nz > 1)
                % �������� �������
                % ���� ����� - �� 2 �����
                % ���� ������� - �� 1 �����
                if rem(Nz, 2) == 0
                    % �����
                    Cz(i) = 0.5 * Nz;
                    Cz(i + 1) = Cz(i) + 1;
                else
                    % �������
                    Cz(i) = round(0.5 * Nz);
                end
            else
                % �� ���� �������, ���� ��������� ����� ����� �����
                % ������� �� {(nz == 1) && (Nz == 1)}
                Cz(i) = 1;
            end
            
        % ������� �����
        elseif (i == nz) && (nz > 1)
            if (nz > 2) && (rem(Nz, 2) == 0)
                Cz(i + 1) = Nz;
            else
                Cz(i) = Nz;
            end
            
        % ���������� �����
        elseif (i == round(0.5 * nz)) && (nz > 2)
            % ���� ��� ������� ����� 2-� �����
            % �������� �������
            % ���� ����� - �� 2 �����
            % ���� ������� - �� 1 �����
            if rem(Nz, 2) == 0
                % �����
                Cz(i) = 0.5 * Nz;
                Cz(i + 1) = Cz(i) + 1;
            else
                % �������
                Cz(i) = round(0.5 * Nz);
            end
        elseif (i > 1) && (i < round(0.5 * nz))
            % ����� �������� �������
            Cz(i) = floor((Nz - 1) * (i - 1) / (nz - 1)) + 1;
            % (Nz - 1) - ������� �������, � �� �������
        elseif (i > round(0.5 * nz)) && (i < nz) && (rem(Nz, 2) == 1)
            % ����� �������� ������� ��� ��������� �������� Nz
            Cz(i) = ceil((Nz - 1) * (i - 1) / (nz - 1)) + 1;
            % (Nz - 1) - ������� �������, � �� �������
        else
            % ����� �������� ������� ��� ������� �������� Nz
            Cz(i + 1) = ceil((Nz - 1) * (i - 1) / (nz - 1)) + 1;
            % (Nz - 1) - ������� �������, � �� �������
        end %
    end
else
    % ������� ���
    for i = 1:Nz
        Cz(i) = i;
    end
end

end
