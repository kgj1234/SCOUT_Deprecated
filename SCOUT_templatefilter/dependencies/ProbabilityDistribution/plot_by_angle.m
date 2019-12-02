function vals=plot_by_angle(xwidth,ywidth,r,theta,initial_angle,modelfun,beta)
radius=xwidth*ywidth./sqrt((xwidth*cos(theta+initial_angle)).^2+(ywidth*sin(theta+initial_angle)).^2);
vals=base_function(r./radius,modelfun,beta);
