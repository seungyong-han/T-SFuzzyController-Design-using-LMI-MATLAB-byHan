function dx=rk5_3(x,u,T)
k1=truck_trailer(x,u)*T;
k2=truck_trailer(x+k1*0.5,u)*T;
k3=truck_trailer(x+k2*0.5,u)*T;
k4=truck_trailer(x+k3,u)*T;
dx=x + ((k1+k4)/6+(k2+k3)/3);
