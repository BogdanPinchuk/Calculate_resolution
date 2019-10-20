function [CWx, CWy] = centerweight(Kx, Ky, MTheta, MPhi, Nx, Ny)
% Функція розрахунку центра мас для пікселів "сусідів" та межі

% Kx, Ky - масив значень нумерацій точок з врахуванням сусідів
% MTheta, MPhi - лінійні координати точок
% Nx, Ny - кількість точок які можливо відобразити

% Розрахуємо граничні значення враховуючи точки сусіди

% Визначаємо min і max границь в циклі
for i = 1:Nx
    for j = 1:Ny
        % Перевіряємо чи піксель без сусідів
        if ((Kx(1, i) - Kx(2, i)) == 0) && ((Ky(1, j) - Ky(2, j)) == 0)
            % Записуємо центри без змін
%             CWx{1}(i, j) = MTheta{1}(i, j);
%             CWy{1}(i, j) = MPhi{1}(i, j);
            
            % Записуємо крайні точки
            for k = 1:5
                CWx{k}(i, j) = MTheta{k}(i, j);
                CWy{k}(i, j) = MPhi{k}(i, j);
            end
        else
            % Записуємо крайні точки
            CWx{2}(i, j) = MTheta{2}(Kx(2, i), Ky(1, j));
            CWy{2}(i, j) = MPhi{2}(Kx(2, i), Ky(1, j));
            
            CWx{3}(i, j) = MTheta{3}(Kx(2, i), Ky(2, j));
            CWy{3}(i, j) = MPhi{3}(Kx(2, i), Ky(2, j));
            
            CWx{4}(i, j) = MTheta{4}(Kx(1, i), Ky(2, j));
            CWy{4}(i, j) = MPhi{4}(Kx(1, i), Ky(2, j));
            
            CWx{5}(i, j) = MTheta{5}(Kx(1, i), Ky(1, j));
            CWy{5}(i, j) = MPhi{5}(Kx(1, i), Ky(1, j));
            
            % Коли в пікселя є сусіди
            [CWx{1}(i, j), CWy{1}(i, j)] = peretintoch(CWx{2}(i, j),...
                CWy{2}(i, j), CWx{4}(i, j), CWy{4}(i, j), CWx{5}(i, j),...
                CWy{5}(i, j), CWx{3}(i, j), CWy{3}(i, j));
        end
    end
end

end
