function [LRatheta, LRaphi] = vidsmiztoch(dovzin, Wx, Wy, Nx, Ny)
% Функція розрахунку роздільної РЗ здатністі враховуючи кривизну,
% але апроксимована лінійно

% LRatheta - РЗ вздовж колонок (вздовж напрямку польоту)
% LRaphi - РЗ вздовж рядків (поперек напрямку польоту)
% dovzin - відстань від КА до кожної точки яка описує РЗ пікселя
% Wx, Wy - кутові величини пікселів (повністю)

% Маємо трикутник в якому розраховуєм за теоремою косинусів РЗ
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
