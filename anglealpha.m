function alpha = anglealpha(theta, phi)
% Функція кута відхилення від надиру

% alpha - кут відхилення точки від осі надиру
% theta - проекція кута відхилення на площину вздовж польоту
% phi - проекція кута відхилення на площину поперек польоту
% Всі кути в радіанах

% вибити помилку через нескінченність при куті в 90°
% alpha = atan(sqrt(tan(theta) .^ 2 + tan(phi) .^ 2));

% Додаткові варіанти розрахунку
% theta2 = atan(tan(phi) .* cos(theta));
% phi2 = atan(tan(theta) .* cos(phi));

% alpha = acos(cos(theta) .* cos(theta2));
% alpha = acos(cos(phi) .* cos(phi2));
% alpha = asin(sqrt(sin(theta2) .^ 2 + sin(phi2) .^ 2));

% Тому уточнимо алгоритм
if abs(phi) == (pi / 2)
    % Якщо кут phi = 90°, то
    phi2 = atan(tan(theta) .* cos(phi));
    alpha = acos(cos(phi) .* cos(phi2));
elseif abs(theta) == (pi / 2)
    % Якщо кут theta = 90°, то
    theta2 = atan(tan(phi) .* cos(theta));
    alpha = acos(cos(theta) .* cos(theta2));
elseif (abs(phi) ==(pi / 2)) && (abs(theta) == (pi / 2))
    alpha = pi / 2;
else
    alpha = atan(sqrt(tan(theta) .^ 2 + tan(phi) .^ 2));
end

end
