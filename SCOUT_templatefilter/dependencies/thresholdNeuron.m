function neuron=thresholdNeuron(neuron,prctile_thresh);
A=neuron.A;
data_shape=neuron.imageSize;
for i=1:size(neuron.A,2)
    A1=A(:,i);
    thresh=max(A1)*prctile_thresh;
    A1(A1<thresh)=0;
    
    A1=reshape(A1,data_shape(1),data_shape(2),[]);
    [centroid,~]=calculateCentroid_and_Covariance(A1,data_shape(1),data_shape(2));
    bw=A1>0;
    stats=regionprops(bw,'MinorAxisLength');
    stats={stats.MinorAxisLength};
    stats=cell2mat(stats);
    len=max(stats);
    circle=createCirclesMask(data_shape,centroid,len/1.75);
    A1(circle==0)=0;
    A1=reshape(A1,data_shape(1)*data_shape(2),[]);
    
    A(:,i)=A1;
end
neuron.A=A;