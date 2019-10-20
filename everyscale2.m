function [MTheta, MPhi] = everyscale2(Nx, Ny, MinT, MinP, Mtheta, Mphi,...
    CWx, CWy)
% ‘ункц≥€ розрахунку центра мас дл€ п≥ксел≥в "сус≥д≥в" та меж≥

% MTheta, MPhi - промасштабованы координати точок п≥ксел≥в
% Mtheta, Mphi - координати точок п≥ксел≥в
% MinT, MinP - масштабн≥ коеф≥ц≥Їнти
% Nx, Ny - к≥льк≥сть точок €к≥ можливо в≥добразити
% CWx, CWy - центра мас дл€ "сус≥д≥в" п≥ксел≥в та меж≥

% ¬изначаЇмо необх≥дний масштаб
if MinT == 1
    % ¬изначаЇмо масштаб по ос≥ OY
    % якщо MinP = NaN (невизначен≥сть)
    if isnan(MinP)
        MinSc = 1.0;
    else
        MinSc = MinP;
    end
else
    % ¬изначаЇмо масштаб по ос≥ OX (де MinP == 1)
    % якщо MinT = NaN (невизначен≥сть)
    if isnan(MinT)
        MinSc = 1.0;
    else
        MinSc = MinT;
    end
end

% ћасштабуЇмо
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
