function Divergences=KLDiv_full(neuron1,neuron2,max_dist_construct)
%Creates JS divergence matrix for footprints in neuron1, neuron2. Assigns 1
%if distance between centroids is larger than max_dist_construct

A1=neuron1.A;
A2=neuron2.A;
A1=A1./sum(A1,1);
A2=A2./sum(A2,1);
dist=pdist2(neuron1.centroid,neuron2.centroid);

%A1=reshape(A1,neuron1.imageSize(1),neuron1.imageSize(2),[]);
%A2=reshape(A2,neuron2.imageSize(1),neuron2.imageSize(2),[]);
Divergences=ones(1,size(A1,2)*size(A2,2));
avail_ind=zeros(1,size(A1,2)*size(A2,2));
parfor a=1:length(avail_ind)
    [i,j]=ind2sub([size(A1,2),size(A2,2)],a);
    if dist(i,j)<max_dist_construct
        avail_ind(a)=1;
    end
end
avail_ind=find(avail_ind); 
val=ones(1,length(avail_ind));
parfor a=1:length(avail_ind)
    [i,j]=ind2sub([size(A1,2),size(A2,2)],avail_ind(a));
      
    val(a)=JSDiv(A1(:,i)',A2(:,j)');
   
end
Divergences(avail_ind)=val;
Divergences=reshape(Divergences,size(A1,2),size(A2,2));
Divergences(isnan(Divergences))=max(reshape(Divergences,1,[]),[],'omitnan');