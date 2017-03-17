function y = gen_shape(x,num)
% Generate wave shapes

x = x - floor(x);
y = zeros(size(x));
switch num
    case 1
        loc = find(x<0.5);
        y(loc) = sqrt(0.25^2-(x(loc)-0.25).^2);
        loc = find(x>=0.5);
        y(loc) = -sqrt(0.25^2-(x(loc)-0.75).^2);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 2
        loc = find(x<=0.25);
        y(loc) = x(loc);
        loc = find(x>0.25 & x<=0.75);
        y(loc) = 0.5-x(loc);
        loc = find(x>0.75);
        y(loc) = x(loc)-1;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 3
        loc = find(x<=0.5);
        y(loc) = 1;
        loc = find(x>0.5);
        y(loc) = -1;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 4
        loc = find(x<=0.25);
        y(loc) = 1;
        loc = find(x>0.25 & x<=0.5);
        y(loc) = 3-8*x(loc);
        loc = find(x>0.5 & x<= 0.75);
        y(loc) = -1;
        loc = find(x>0.75);
        y(loc) = 8*x(loc)-7;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 5
        loc = find(x<=0.45);
        y(loc) = x(loc);
        loc = find(x>0.45 & x<=0.55);
        y(loc) = 4.5-9*x(loc);
        loc = find(x>0.55);
        y(loc) = x(loc)-1;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 6
        x = 2*pi*x;
        y = (cos(0.8*cos(x))-tan(x).*sin(0.8*cos(x))).*cos(x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 7
        loc = find(x<=0.3);
        y(loc) = 0.5;
        loc = find(x>0.3 &x<=0.6);
        y(loc) = -0.2;
        loc = find(x>0.6);
        y(loc) = 0.2;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 8
        x = 2*pi*x;
        y = cos(x) + 2*cos(2*x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 9
        x = 2*pi*x;
        y = cos(x) + 3+cos(3*x)+ 2*cos(4*x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 10
        x = 2*pi*x;
        y = cos(x) + 2*cos(3*x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    case 11
        x = 2*pi*x;
        y = cos(x) + 2*cos(4*x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
end