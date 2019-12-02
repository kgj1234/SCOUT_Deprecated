function centroid=calculateCentroid(A_individual,height,width)
A_individual=reshape(A_individual,[height,width]);
centroid=zeros(1,2);
for i=1:height
    for j=1:width
        centroid=centroid+[j,i]*A_individual(i,j);
    end
end
centroid=centroid/sum(sum(A_individual));

