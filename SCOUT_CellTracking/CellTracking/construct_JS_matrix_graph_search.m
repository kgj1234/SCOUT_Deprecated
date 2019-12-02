function JS_matrices=construct_JS_matrix_graph_search(neurons,max_dist_construct)

%Construct JS distance between all neuron pairs extracted from
%different recordings
%inputs
%   neurons: cell array of extracted neural activity from each recording,
%       Sources2D

%outputs
%   JS_matrices: cell array of JS distance matrices
%Author: Kevin Johnston, University of California, Irvine


total=0;
clear display_progress_bar
display_progress_bar('Computing JS matrices: ',false);

for i=1:length(neurons)-1
    for j=i:length(neurons)
        JS_matrices{i,j}=KLDiv_full(neurons{i},neurons{j},max_dist_construct);
        
        total=total+1;
        display_progress_bar((total/((length(neurons)+1)*(length(neurons))/2-1)*100),false);
    end
end
display_progress_bar(' Completed',false);
display_progress_bar('terminated',true);