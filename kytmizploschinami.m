function deltaAlpha = kytmizploschinami(fi, theta, phi, H, Rk)
% Функція визначення кута між двома площинами, 
% по прямій яка задана кутом від напрямку польоту

% fi - кут між напрямком польоту та лінією по якій необхідно знайти кут між двома площинами
% theta - Кут візування в напрямку польоту
% phi - Кут візування поперек польоту

% Визначаємо кут відхилення КА
alphaKA = anglealpha(theta, phi);

% Визначаємо відповідний йому кут в центрі Землі (кут між точкою надиру і 
% точкою відхилення КА на поверхні Землі)
thetaZ = angleZemli(H, Rk, alphaKA);

% Подібно до полярних координат, визначаємо кут між напрямком перпендикулярно польоту і проекційною
% лінією (від точки нахилу КА на Землі до точки надиру) лінії нахилу КА
if ((theta == 0) && ((0 <= phi) || (phi < 0.5 * pi)))
    phi_1 = 0;
elseif ((theta == 0) && ((-0.5 * pi < phi) || (phi < 0)))
    phi_1 = pi;
elseif ((theta ~= 0) && (abs(theta) < 0.5 * pi) && (phi == 0))
    phi_1 = sign(theta) * 0.5 * pi;
elseif (((0 < abs(theta)) || (abs(theta) < 0.5 * pi)) &&...
        ((0 < phi) || (phi < 0.5 * pi)))
    phi_1 = atan(tan(theta) / tan(phi));
elseif (((0 < abs(theta)) || (abs(theta) < 0.5 * pi)) &&...
        ((-0.5 * pi < phi) || (phi < 0)))
    phi_1 = pi + atan(tan(theta) / tan(phi));
else
    phi_1 = 0;
end

% Кут між двома площинами
deltaAlpha = 2 * asin(sin(phi_1 + fi) * sin(0.5 * thetaZ));

end