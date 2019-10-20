function [Theta, Phi] = naklonLok2(theta, phi, Wx, Wy, ind)
% Функція розрахунку нахилу КА

% Theta - theta1
% Phi - phi1
% theta, phi - стандарт
% Wx - theta0
% Wy - phi0
% ind - індекс який визначає як рухатиметься система
% якщо 0 - то спочатку по тангажу, а тоді по крену
% якщо 1 - то спочатку по крену, а тоді по тангажу

% Попередній розрахунок
% theta_2 = atan(tan(phi) .* cos(theta));
% phi_2 = atan(tan(theta) .* cos(phi));

Wx_2 = atan(tan(Wy) .* cos(Wx));
Wy_2 = atan(tan(Wx) .* cos(Wy));

alpha = atan(sqrt(tan(theta) .^ 2.0 + tan(phi) .^ 2.0));
Walpha = atan(sqrt(tan(Wx) .^ 2.0 + tan(Wy) .^ 2.0));

alpha2 = acos(cos(alpha) .* cos(Walpha) .* (1.0 - tan(theta) .* tan(Wx)...
        - tan(phi) .* tan(Wy)));

if ind == 0
   % якщо ind = 0
%    dphi = atan((tan(Wy) .* cos(Wx)) ./ cos(Wx + phi_2));
%    Phi = phi + dphi;
%    Theta = atan((tan(Wx + phi_2) .* cos(dphi)) ./ cos(Phi));
    alpha1 = acos(cos(theta) .* cos(Walpha) .* (1.0 - tan(theta) .* tan(Wx)));
    
    theta1 = acos(cos(alpha1) ./ cos(Wx_2));
    phi1 = atan(tan(Wx_2) ./ cos(theta1));
    
    phi2 = acos(cos(alpha2) .* cos(phi1) ./ cos(alpha1));
%     phi2 = acos(cos(alpha2) ./ cos(atan(tan(theta1) .* cos(phi1))));
    theta2 = atan(cos(alpha1) .* tan(theta1) ./ cos(alpha2));
%     theta2 = atan(tan(theta1) .* cos(phi1) ./ cos(phi2));
else
   % якщо ind = 1
%    dtheta = atan((tan(Wx) .* cos(Wy)) ./ cos(Wy + theta_2));
%    Theta = theta + dtheta;
%    Phi = atan((tan(Wy + theta_2) .* cos(dtheta)) ./ cos(Theta));
    alpha1 = acos(cos(phi) .* cos(Walpha) .* (1.0 - tan(phi) .* tan(Wy)));
    
    phi1 = acos(cos(alpha1) ./ cos(Wy_2));
    theta1 = atan(tan(Wy_2) ./ cos(phi1));
    
    theta2 = acos(cos(alpha2) .* cos(theta1) ./ cos(alpha1));
    phi2 = atan(cos(alpha1) .* tan(phi1) ./ cos(alpha2));
end

Theta = real(theta2);
Phi = real(phi2);

end
