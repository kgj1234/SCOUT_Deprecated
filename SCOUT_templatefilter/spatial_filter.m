function [neuron,KL]=spatial_filter(neuron,KL_thresh,data_shape,constraint_type,trim,gSizMax,gSizMin,filter,update_method,div, threshold_per)
n=size(neuron.C,1);
%constraint_type: 'bound', 'prc'

if ~exist('constraint_type','var')||isempty(constraint_type)
    constraint_type='prc';
end
if ~exist('threshold_per','var')||isempty(threshold_per)
    threshold_per=.1;
end
if ~exist('update_method','var')||isempty(update_method)
    update_method='threshold';
end
if ~exist('div','var')||isempty(div)
    div='JS';
end
if ~exist('gSizMin','var')||isempty(gSizMin);
    gSizMin=7;
end

sigma=(gSizMin-1)/6;
h=fspecial('gaussian',[gSizMin,gSizMin],sigma);


if size(data_shape,1)>1;
    data_shape=squeeze(data_shape');
end
det_cov=[];
eig_ratio=[];
KL=[];
for i=1:size(neuron.A,2);
    A{i}=full(neuron.A(:,i));
end
ind_del=[];
n=length(A);
KL1=2*ones(1,length(A));
KL2=2*ones(size(KL1));
mass1=ones(size(KL1));
mass2=ones(size(KL2));
parfor i=1:n
    
    try
        A1=A{i};
        
        A1=reshape(A1,data_shape);
        
        
        
        curr_thresh=min(threshold_per,.85);
        A3=A1;
        
        while curr_thresh<=.85
            
            %A1=imgaussfilt(A1,1.3);
            thresh_ind=A1<curr_thresh*max(max(A1));
            A1(thresh_ind)=0;
            bw=A1>0;
            stats= regionprops(full(bw),'MajorAxisLength','MinorAxisLength');
            MajorAxisLength={stats.MajorAxisLength};
            MinorAxisLength={stats.MinorAxisLength};
            MajorAxisLength=cell2mat(MajorAxisLength);
            MinorAxisLength=cell2mat(MinorAxisLength);
            MajorAxisLength=max(MajorAxisLength);
            MinorAxisLength=max(MinorAxisLength);
            if MinorAxisLength<gSizMin
                break
            end
            if MajorAxisLength<.95*gSizMax
                break
            end
            
            curr_thresh=curr_thresh+.05;
        end
        
        
        if MinorAxisLength<gSizMin || curr_thresh>=.89
            continue;
        end
        [centroid,cov]=calculateCentroid_and_Covariance(A1,data_shape(1),data_shape(2));
        sum_A=sum(sum(A1));
        A1=A1/sum_A;
 
        
        if trim
        [A1_comp,A2_comp]=construct_comparison_footprint_ellipse(A1,centroid,data_shape);
        %[A1_comp,A2_comp]=construct_comparison_footprint_gaussian(centroid,cov,data_shape);
        
        KL1(i)=JSDiv(reshape(A1,1,[]),reshape(A1_comp,1,[]));
        KL2(i)=JSDiv(reshape(A1,1,[]),reshape(A2_comp,1,[]));
        else
            [A1_comp,A2_comp]=construct_comparison_footprint_ellipse(A3,centroid,data_shape);
        %[A1_comp,A2_comp]=construct_comparison_footprint_gaussian(centroid,cov,data_shape);
        
         KL1(i)=JSDiv(reshape(A3,1,[]),reshape(A1_comp,1,[]));
        KL2(i)=JSDiv(reshape(A3,1,[]),reshape(A2_comp,1,[]));
        end
        best_choice=1;
        if KL2(i)<KL1(i)
            best_choice=2;
        end
        if best_choice==1
            if isequal(update_method,'threshold')
                thresh_indices=A1_comp<.1*max(max(A1_comp));
                A3(thresh_indices)=0;
            end
            %Bayesian Updata
            %A1=A1.*A1_comp;
            %A1=A1/sum(sum(A1));
            %A1=A1*sum_A;
            
            %A1=(A3+A1_comp)/2;
            
        else
            if isequal(update_method,'threshold')
                thresh_indices=A2_comp<.1*max(max(A2_comp));
                A3(thresh_indices)=0;
            end
            %A1=A1.*A2_comp;
            %A1=A1/sum(sum(A1));
            %A1=A1*sum_A;
            
            %A1=(A3+A2_comp)/2;
            
        end
        A1=A1*sum_A;
        A1=imfilter(A1,h);
        
        if ~trim
            A1=A3;
        end
        bw=A1>0;
            stats= regionprops(full(bw),'MajorAxisLength','MinorAxisLength');
            MajorAxisLength={stats.MajorAxisLength};
            MinorAxisLength={stats.MinorAxisLength};
            MajorAxisLength=cell2mat(MajorAxisLength);
            MinorAxisLength=cell2mat(MinorAxisLength);
            MajorAxisLength=max(MajorAxisLength);
            MinorAxisLength=max(MinorAxisLength);
        Axis_Ratio(i)=MajorAxisLength/MinorAxisLength;
        
    catch
        mass1(i)=1;
        mass2(i)=1;
        KL1(i)=KL_thresh+2;
        KL2(i)=KL_thresh+2;
 %       max_diam(i)=0;
 %       min_diam(i)=0;
        
    end
    
    
    
    
    A{i}=reshape(A1,1,[]);
    
end

for i=1:length(A)
    neuron.A(:,i)=reshape(A{i},[],1);
end

KL1=[KL1,zeros(1,size(neuron.A,2)-length(KL1))];

KL2=[KL2,zeros(1,size(neuron.A,2)-length(KL2))];

mass1=[mass1,zeros(1,size(neuron.A,2)-length(mass1))];

mass2=[mass2,zeros(1,size(neuron.A,2)-length(mass2))];


%max_diam=[max_diam,zeros(1,size(neuron.A,2)-length(max_diam))];
%max_diam(ind_del)=[];
%min_diam=[min_diam,zeros(1,size(neuron.A,2)-length(min_diam))];
%min_diam(ind_del)=[];
neuron.delete(ind_del)
neuron.A=sparse(neuron.A);
if filter==true
    KL=min([reshape(KL1,1,[]);reshape(KL2,1,[])]);
    
    
    
    
    
    mass=min([reshape(mass1,1,[]);reshape(mass2,1,[])]);
    
    
    
    
    
    
    if isequal(constraint_type,'prc')
        try
            KL=(KL-mean(KL))/std(KL);
            KL_thresh=norminv(KL_thresh);
        end
    end
    
    indices=find(KL>KL_thresh);
 
else
    indices=[];
end


for i=1:size(neuron.A,2)
    if sum(neuron.A(:,i))<=0
        indices=[indices,i];
    end
end
neuron.delete(indices);
KL(indices)=[];

m=size(neuron.C,1);
disp(['Neurons Deleted by Spatial Filter: ',num2str(n-m)]);
end
