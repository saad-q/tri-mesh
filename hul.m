function [X,Y] = hul(xx,yy)



cent = [xx(end) yy(end)];

L = length(xx) - 1;
Ang = zeros(L,1);

X = Ang; Y = X;


for j = 1:L
    a = atan(abs(yy(j) - cent(2))/abs(xx(j) - cent(1)));
    
    if yy(j) - cent(2) >=0 && xx(j) - cent(1) >=0
        Ang(j) = a;
    elseif yy(j) - cent(2) >=0 && xx(j) - cent(1) < 0
        Ang(j) = pi - a;
    elseif yy(j) - cent(2) < 0 && xx(j) - cent(1) < 0
        Ang(j) = pi + a;
    else Ang(j) = 2*pi-a;
    end
end
        

for i = 1:L
    [~ , k] = min(Ang);
    X(i) = xx(k);
    Y(i) = yy(k);
    
    Ang(k) = 100;
end
    





end




