function [EllipseX,EllipseY]=createEllipse()
    majorAxisX1=2.82;
    majorAxisX2=-2.82;
    majorAxisY1=0.000;
    majorAxisY2=0.00;
    e=0.855;
    a = 1/2*sqrt((majorAxisX2-majorAxisX1)^2+(majorAxisY2-majorAxisY1)^2); 
    b = a*sqrt(1-e^2); t = linspace(0,2*pi);
    X = a*cos(t); Y = b*sin(t);
    w = atan2(majorAxisY2-majorAxisY1,majorAxisX2-majorAxisX1); 
    EllipseX = (majorAxisX1+majorAxisX2)/2 + X*cos(w) - Y*sin(w);
    EllipseY = (majorAxisY1+majorAxisY2)/2 + X*sin(w) + Y*cos(w);
end