function aligned_probabilities=construct_combined_probabilities(aligned_neurons,probabilities,pair_aligned,spat_aligned,min_prob,method,penalty)
%methods: min_intergroup, mean_intergroup,
%min_max,mean_max,total_connected,individual_connected
%default min_max

if ~exist('method','var')||isempty(method)
    method='individual_connected';
end
if ~exist('penalty','var')||isempty(penalty)
    penalty=min_prob/2;
end

aligned_probabilities=zeros(size(aligned_neurons,1),1);


for i=1:size(aligned_neurons,1)
    if sum(~iszero(aligned_neurons(i,:)))>0
        
        
        
        total_prob=0;
        min_curr_prob=1;
        max_group=[];
        individual_connected=[];
        total_exceeding=0;
        total_pairs=0;
        for j=1:size(aligned_neurons,2)
            index_prob=0;
            connected=0;
            
            for k=[1:j-1,j+1:size(aligned_neurons,2)]
                if j>k
                    index1=k;
                    index2=j;
                else
                    index1=j;
                    index2=k;
                end
                
                if ~iszero(aligned_neurons(i,index1))&~iszero(aligned_neurons(i,index2))
                    
                    if index2==index1+1
                        ind=find(pair_aligned{index1}(:,1)==aligned_neurons(i,index1)&pair_aligned{index1}(:,2)==aligned_neurons(i,index2));
                        if isempty(ind)
                            total_prob=total_prob+penalty;
                        else
                            total_prob=total_prob+probabilities{index1,index2}(ind);
                            index_prob=max(index_prob,probabilities{index1,index2}(ind));
                            total_exceeeding=total_exceeding+1;
                            connected=connected+1;
                            if probabilities{index1,index2}(ind)<min_curr_prob
                                min_curr_prob=probabilities{index1,index2}(ind);
                            end
                        end
                        
                        total_pairs=total_pairs+1;
                        
                        
                    else
                          
                        ind=find(spat_aligned{index1,index2}(:,1)==aligned_neurons(i,index1)&spat_aligned{index1,index2}(:,2)==aligned_neurons(i,index2));
                        if isempty(ind)
                            total_prob=total_prob+penalty;
                        else
                            total_prob=total_prob+probabilities{index1,index2}(ind);
                            index_prob=max(index_prob,probabilities{index1,index2}(ind));
                            total_exceeding=total_exceeding+1;
                            connected=connected+1;
                            if probabilities{index1,index2}(ind)<min_curr_prob
                                min_curr_prob=probabilities{index1,index2}(ind);
                            end
                        end
                        total_pairs=total_pairs+1;
                        
                        
                    end
                end
            end
        max_group=[max_group,index_prob];    
        individual_connected=[individual_connected,connected];
        end
        max_group=max_group(~iszero(aligned_neurons(i,:)));
        if isequal(method,'min_max');
            aligned_probabilities(i)=min(max_group);
        elseif isequal(method,'mean_max');
            aligned_probabilities(i)=mean(max_group);
        elseif isequal(method,'mean_intergroup')
            aligned_probabilities(i)=total_prob/total_pairs;
        elseif isequal(method,'total_connected')
            aligned_probabilities(i)=total_exceeding/total_pairs;
        elseif isequal(method,'individual_connected')
            if sum(~iszero(aligned_neurons(i,:)))>1&sum(~iszero(individual_connected))>0
                %Adding the additional portion allows for greater
                %differentiation when multiple similar probabilities are
                %obtained.
                ind=find(~iszero(aligned_neurons(i,:)));
                aligned_probabilities(i)=min(individual_connected(ind))/(sum(~iszero(aligned_neurons(i,:)))-1)+total_prob/total_pairs*10^(-2);
            elseif sum(~iszero(aligned_neurons(i,:)))==1
                aligned_probabilities(i)=1;
            else
                aligned_probaiblities(i)=0;
            end
        else
            aligned_probabilities(i)=min_curr_prob;
        end
        
        
    else
        aligned_probabilities(i)=0;
    end
end