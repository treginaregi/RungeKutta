%Método runge kutta de cuarto orden
function [y] = RungeKutta(f,t0,tf,h,yn)
    n=(tf-t0)/h;
    y(1)=yn;
    for i=1:n
        tn=t0+i*h;
        k1=h*f(tn,y(i));
        k2=h*f(tn+0.5*h,y(i)+0.5*k1);
        k3=h*f(tn+0.5*h,y(i)+0.5*k2);
        k4=h*f(tn+h,y(i)+k3);
        y(i+1)=y(i)+(1/6)*(k1+2*k2+2*k3+k4);
    end
end