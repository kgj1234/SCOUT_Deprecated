function [A1_comp,A2_comp]=construct_comparison_footprint_ellipse(A1,centroid,data_shape)
    A1=reshape(A1,data_shape);
    bw=A1>0;
    CC=bwconncomp(bw);
    stats=CC.PixelIdxList;
    curr_rem_ind=[];
    num_elements=[];
    for l=1:length(stats)
        num_elements(l)=length(stats{l});
    end
    for l=1:length(stats)
        if length(stats{l})~=max(num_elements)
            curr_rem_ind=[curr_rem_ind;stats{l}];
        end
    end
    A2=A1;
    A2(curr_rem_ind)=0;
    bw=A2>0;
    stats= regionprops(full(bw),'MajorAxisLength','MinorAxisLength','Orientation');
    MajorAxisLength=stats.MajorAxisLength;
    
    MinorAxisLength=stats.MinorAxisLength;
   
    Orientation=stats.Orientation;
   
    %Estimate appropriate power to use
    x=round(MajorAxisLength*cos(Orientation*pi/180));
    y=round(MajorAxisLength*sin(Orientation*pi/180));
    edge_vec=centroid+[x,y];
    x_coord=[centroid(1),edge_vec(1)];
    y_coord=[centroid(2),edge_vec(2)];
    c = improfile(A1,x_coord,y_coord);
    
    c=c/max(c);
    c(c==0)=[];
    c(isnan(c))=[];
    c=c';
    c=[c,0];
    x=linspace(0,1,length(c));
    error=[nan,nan];
    x=x(1:end-1);
    c=c(1:end-1);
    warning('off','all')
%     try
%         modelfun{1}=@(b,x) exp(-b(1)*(x+b(2)))./(1+exp(-b(1)*(x+b(2))))+b(3);
%         beta{1}=nlinfit(x,c,modelfun{1},[1,0,.1]);
%         vals=modelfun{1}(beta{1},x);
%         error(1)=sqrt(sum((vals-c).^2)/length(vals));
%     
%     end
%     try
        modelfun{2}=@(b,x) (1-x).^(1/b(1));
        beta{2}=nlinfit(x,c,modelfun{2},[2.2]);
        vals=modelfun{2}(beta{2},x);
        %error(2)=sqrt(sum((vals-c).^2)/length(vals));
    %end
    warning('on','all')
%[~,I]=min(error);
I=2;
beta=beta{I};
modelfun=modelfun{I};
  
        
    
    
    A1_comp=plot_ellipse((MajorAxisLength/2),(MinorAxisLength/2),modelfun,beta,Orientation*pi/180,centroid,data_shape);
    A2_comp=plot_ellipse((MajorAxisLength/2),(MinorAxisLength/2),modelfun,beta,Orientation*pi/180+pi/2,centroid,data_shape);
    