function A=plot_ellipse(xwidth,ywidth,modelfun,beta,angle,centroid,data_shape)
[t,r] = meshgrid(linspace(0,2*pi,361),linspace(0,max(xwidth,ywidth),ceil(2*max(xwidth,ywidth))));
%t=reshape(t,1,[]);
%r=reshape(r,1,[]);
Z=plot_by_angle(xwidth,ywidth,r,t,angle,modelfun,beta);
Z(abs(imag(Z))>0)=0;


[X,Y]=pol2cart(t,r);
X=round(X);
Y=round(Y);




X=reshape(X,1,[])+round(centroid(1));
Y=reshape(Y,1,[])+round(centroid(2));
X(X>data_shape(2))=data_shape(2);
X(X<1)=1;
Y(Y>data_shape(1))=data_shape(1);
Y(Y<1)=1;
A=zeros(data_shape);
ind=sub2ind(data_shape,Y,X);
for i=1:length(ind)
    if ind(i)<=0
        ind=sub2ind(data_shape,round(centroid(1)),round(centroid(2)));
    end
        
    same_index=ind==ind(i);
    Z(same_index)=max(Z(same_index));
end


A(ind)=Z;
A=A/sum(sum(A));


