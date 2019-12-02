function [neuron,Path_Storage,neurons,links]=cellTracking_SCOUT(neurons,varargin)
%Author: Kevin Johnston
%Main function for cell registration between sessions
%Inputs
%neurons (cell of class Sources2D objects containing neural data from each
%           session, A,C,imageSize required entries)
%Optional Inputs
%  links: extractions of connecting sesssions, required for correlation
%               registration. Cell of Sources2D structs
%  overlap: temporal overlap between connecting sessions
%  max_dist: max distance between registered neurons between sessions
%  weights: weights for linking methods. 4 element numeric vector. Order correlation, JS, overlap, centroid dist
%  chain_prob: Total probability threshold for neuron chains, numeric value between 0 and 1
%  corr_thresh Probability threshold for correlation link
%               between neurons, numeric value between 0 and 1. Usually 0
%               unless binary_corr is true
%  patch: Patch size if neuron density varies significantly across the
%               recording
%  register_sessions: boolean indicating whether to register sessions.
%       Default true
%  registration_type: 'align_to_base' or 'consecutive', Determines if
%       sessions are registered to base session, or registered
%       consecutively.
%  registration_method: 'affine' or 'non-rigid', method for session
%       registration
%  registration_template: Use spatial positions of neurons 'spatial', or
%       correlation map 'correlation' for registration
%  use_corr: boolean indicating whether to use correlation cell registration
%           requires links
%  use_spat: boolean indicating whether to use spatial cell registration.
%       This should be false only if resources are low, or registration
%       cannot be guaranteed between all sessions
%   max_gap: Integer indicating largest gap between sessions for cell
%       registration. 0 indicates cells must appear in each session
%   probability_assignment_method: either 'percentile', 'gmm',
%       'gemm','glmm' 'default' or 'Kmeans'
%   base: base session for alignment
%   min_prob: minimum probability for identification between
%   sessions
%   binary_corr: treat correlation on connecting recordings as binary
%   variable
% Outputs
%  neuron: Sources2D containing neural activity extracted through the set
%           of sessions
%  Path_Storage: Registration Indices Per Session for each identified cell
%%
%Assign variables
p = inputParser;
p.CaseSensitive=false;
defaultWeights = [1,5,5,5,1,1];
defaultLinks=[];
defaultProbThresh=.5;
defaultCorrThresh=.7;
defaultSinglecorr=false;
defaultPatch=[];
defaultRegistration='align_to_base';
defaultRegistrationtemplate='spatial';
defaultRegistrationmethod='affine';
defaultCorr=true;
defaultOverlap=0;
defaultMinProb=.5;
defaultSpat=true;
defaultBincor=false;
defaultMaxgap=0;
defaultMaxdist=15;
defaultRegister=true;
defaultProbmethod='default';
defaultBase=ceil(length(neurons)/2);
addRequired(p,'neurons');
addOptional(p,'links',defaultLinks);
addOptional(p,'overlap',defaultOverlap);
addParameter(p,'weights',defaultWeights);
addParameter(p,'single_corr',defaultSinglecorr);
addParameter(p,'chain_prob',defaultProbThresh,@isnumeric);
addParameter(p,'corr_thresh',defaultCorrThresh,@isnumeric);
addParameter(p,'registration_type',defaultRegistration,@isstr);
addParameter(p,'min_prob',defaultMinProb,@isnumeric);
addParameter(p,'patch',defaultPatch);
addParameter(p,'use_corr',defaultCorr);
addParameter(p,'use_spat',defaultSpat);
addParameter(p,'register_sessions',defaultRegister);
addParameter(p,'registration_template',defaultRegistrationtemplate,@isstr);
addParameter(p,'registration_method',defaultRegistrationmethod,@isstr);
addParameter(p,'max_dist',defaultMaxdist,@isnumeric);
addParameter(p,'max_gap',defaultMaxgap,@isnumeric);
addParameter(p,'probability_assignment_method',defaultProbmethod,@isstr);
addParameter(p,'base',defaultBase,@isnumeric);
addParameter(p,'binary_corr',defaultBincor);
parse(p,neurons,varargin{:});
for i=1:length(p.Parameters)
    val=getfield(p.Results,p.Parameters{i});
    eval([p.Parameters{i},'=val']);
end
if isempty(links)
    weights(1)=0;
    overlap=0;
    corr_thresh=[];
end
if isempty(neurons{1}.Cn)
    registration_template='spatial';
end
if weights(1)==0
    links=[];
    overlap=0;
    corr_thresh=[];
end
%% Copy neurons and links to new memory locations, standardize FOV size for each session


%Sources2D is mutable, copy a new version, copy a version without trimmed
%neurons
neurons1=neurons;
links1=links;
for i=1:length(neurons)
    neurons1{i}=neurons{i}.copy();
end
clear neurons
neurons=neurons1;
clear neurons1;
for i=1:length(neurons)
    neurons1{i}=neurons{i}.copy();
    neurons{i}=thresholdNeuron(neurons{i},.4);
end
for i=1:length(links)
    links1{i}=links1{i}.copy();
end
clear links
links=links1;
clear links1;
for i=1:length(links)
    links1{i}=links{i}.copy();
    links{i}=thresholdNeuron(links{i},.4);
end
%Normalize FOV for each session
for i=1:length(neurons)
    
    
    max_dims1(i)=neurons{i}.imageSize(1);
    max_dims2(i)=neurons{i}.imageSize(2);
end
for i=1:length(links)
    max_dims1(end+1)=links{i}.imageSize(1);
    max_dims2(end+1)=links{i}.imageSize(2);
end
max_dims1=max(max_dims1);
max_dims2=max(max_dims2);
for i=1:length(neurons)
    curr_dim1=max_dims1-neurons{i}.imageSize(1);
    curr_dim2=max_dims2-neurons{i}.imageSize(2);
    neurons{i}.A=reshape(neurons{i}.A,neurons{i}.imageSize(1),neurons{i}.imageSize(2),[]);
    neurons{i}.A=[neurons{i}.A;zeros(curr_dim1,size(neurons{i}.A,2),size(neurons{i}.A,3))];
    neurons{i}.A=[neurons{i}.A,zeros(size(neurons{i}.A,1),curr_dim2,size(neurons{i}.A,3))];
    neurons{i}.imageSize=[max_dims1,max_dims2];
    neurons{i}.A=reshape(neurons{i}.A,max_dims1*max_dims2,[]);
    try
        neurons{i}.Cn=[neurons{i}.Cn;zeros(curr_dim1,size(neurons{i}.Cn,2))];
        neurons{i}.Cn=[neurons{i}.Cn,zeros(size(neurons{i}.Cn,1),curr_dim2)];
    catch
    end
end
if ~isempty(links)
    for i=1:length(links)
        curr_dim1=max_dims1-links{i}.imageSize(1);
        curr_dim2=max_dims2-links{i}.imageSize(2);
        links{i}.A=reshape(links{i}.A,links{i}.imageSize(1),links{i}.imageSize(2),[]);
        links{i}.A=[links{i}.A;zeros(curr_dim1,size(links{i}.A,2),size(links{i}.A,3))];
        links{i}.A=[links{i}.A,zeros(size(links{i}.A,1),curr_dim2,size(links{i}.A,3))];
        links{i}.imageSize=[max_dims1,max_dims2];
        links{i}.A=reshape(links{i}.A,max_dims1*max_dims2,[]);
        try
            links{i}.Cn=[links{i}.Cn;zeros(curr_dim1,size(links{i}.Cn,2))];
            links{i}.Cn=[links{i}.Cn,zeros(size(links{i}.Cn,1),curr_dim2)];
        catch
        end
    end
end


%% Register Sessions


if register_sessions
    
    clear display_progress_bar
    display_progress_bar('Aligning Recordings: ',false);
    display_progress_bar(0,false);
    if isequal(registration_template,'correlation')
        base_template=neurons{base}.Cn;
        base_template(base_template<.9)=0;
    else
        base_template=max(reshape(neurons{base}.A./max(neurons{base}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
    end
    templates{base}=base_template;
    if isequal(registration_type,'align_to_base')
        for i=[1:base-1,base+1:length(neurons)]
            if isequal(registration_template,'correlation')
                template=neurons{i}.Cn;
                template(template<.9)=0;
            else
                template=max(reshape(neurons{i}.A./max(neurons{i}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
            end
            
            registration1=registration2d(base_template,template,'transformationModel','translation');
            template=imtranslate(template,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            registration2=registration2d(base_template,template,'transformationModel',registration_method);
            templates{i}=deformation(template,registration2.displacementField,registration2.interpolation);
            temp_A=reshape(neurons{i}.A,neurons{i}.imageSize(1),neurons{i}.imageSize(2),[]);
            parfor j=1:size(neurons{i}.A,2)
                temp_A(:,:,j)=imtranslate(temp_A(:,:,j),-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                temp_A(:,:,j)=deformation(temp_A(:,:,j),registration2.displacementField,registration2.interpolation);
            end
            neurons{i}.A=reshape(temp_A,neurons{i}.imageSize(1)*neurons{i}.imageSize(2),[]);
            neurons{i}.A(neurons{i}.A<10^(-6))=0;
            if ~isempty(neurons{i}.Cn)
            neurons{i}.Cn=imtranslate(neurons{i}.Cn,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            neurons{i}.Cn=deformation(neurons{i}.Cn,registration2.displacementField,registration2.interpolation);
            end
            %    neurons{i}=thresholdNeuron(neurons{i},.3);
            display_progress_bar(i/length(neurons)*100,false);
        end
    else
        total=1;
        templates{base}=base_template;
        for i=base+1:length(neurons)
            if isequal(registration_template,'correlation')
                template_curr=neurons{i}.Cn;
                template_curr(template_curr<.9)=0;
                template_prev=neurons{i-1}.Cn;
                template_prev(template_prev<.9)=0;
            else
                template_curr=max(reshape(neurons{i}.A./max(neurons{i}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
                template_prev=max(reshape(neurons{i-1}.A,neurons{1}./max(neurons{i-1}.A,[],1).imageSize(1),neurons{1}.imageSize(2),[]),[],3);
            end
            
            registration1=registration2d(template_prev,template_curr,'transformationModel','translation');
            template=imtranslate(template_curr,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            registration2=registration2d(template_prev,template,'transformationModel',registration_method);
            templates{i}=deformation(template,registration2.displacementField,registration2.interpolation);
            temp_A=reshape(neurons{i}.A,neurons{i}.imageSize(1),neurons{i}.imageSize(2),[]);
            parfor j=1:size(neurons{i}.A,2)
                temp_A(:,:,j)=imtranslate(temp_A(:,:,j),-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                temp_A(:,:,j)=deformation(temp_A(:,:,j),registration2.displacementField...
                    ,registration2.interpolation);
            end
            neurons{i}.A=reshape(temp_A,neurons{i}.imageSize(1)*neurons{i}.imageSize(2),[]);
            %non-rigid registration can cause spurious non zero values
            neurons{i}.A(neurons{i}.A<10^(-6))=0;
            if ~isempty(neurons{i}.Cn)
            neurons{i}.Cn=imtranslate(neurons{i}.Cn,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            neurons{i}.Cn=deformation(neurons{i}.Cn,registration2.displacementField,registration2.interpolation);
            end
            total=total+1;
            display_progress_bar(total/length(neurons)*100,false);
        end
        for i=base-1:-1:1
            if isequal(registration_template,'correlation')
                template_curr=neurons{i}.Cn;
                template_curr(template_curr<.9)=0;
                template_prev=neurons{i+1}.Cn;
                template_prev(template_prev<.9)=0;
            else
                template_curr=max(reshape(neurons{i}.A./max(neurons{i}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
                template_prev=max(reshape(neurons{i+1}.A./max(neurons{i+1}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
            end
            
            registration1=registration2d(template_prev,template_curr,'transformationModel','translation');
            template=imtranslate(template_curr,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            registration2=registration2d(template_prev,template,'transformationModel',registration_method);
            templates{i}=deformation(template,registration2.displacementField,registration2.interpolation);
            temp_A=reshape(neurons{i}.A,neurons{i}.imageSize(1),neurons{i}.imageSize(2),[]);
            parfor j=1:size(neurons{i}.A,2)
                temp_A(:,:,j)=imtranslate(temp_A(:,:,j),-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                temp_A(:,:,j)=deformation(temp_A(:,:,j),registration2.displacementField...
                    ,registration2.interpolation);
            end
            neurons{i}.A=reshape(temp_A,neurons{i}.imageSize(1)*neurons{i}.imageSize(2),[]);
            %non-rigid registration can cause spurious non zero values
            neurons{i}.A(neurons{i}.A<10^(-6))=0;
            if ~isempty(neurons{i}.Cn)
            neurons{i}.Cn=imtranslate(neurons{i}.Cn,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
            neurons{i}.Cn=deformation(neurons{i}.Cn,registration2.displacementField,registration2.interpolation);
            end
            total=total+1;
            display_progress_bar(total/length(neurons)*100,false);
        end
    end
    
    display_progress_bar(' Completed',false);
    display_progress_bar('terminating',true)
    
    
    
    clear display_progress_bar
    if ~isempty(links)
        display_progress_bar('Aligning Connecting Recordings: ',false);
        display_progress_bar(0,false);
        if isequal(registration_type,'align_to_base')
            
            for i=1:length(links)
                if isequal(registration_template,'correlation')
                    template=links{i}.Cn;
                    template(template<.9)=0;
                else
                    template=max(reshape(links{i}.A./max(links{i}.A,[],1),neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]),[],3);
                end
                
                
                registration1=registration2d(base_template,template,'transformationModel','translation');
                template=imtranslate(template,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                registration2=registration2d(base_template,template,'transformationModel',registration_method);
                template_links{i}=deformation(template,registration2.displacementField,registration2.interpolation);
                temp_A=reshape(links{i}.A,neurons{1}.imageSize(1),neurons{1}.imageSize(2),[]);
                parfor j=1:size(links{i}.A,3)
                    temp_A(:,:,j)=imtranslate(temp_A(:,:,j),-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                    temp_A(:,:,j)=deformation(temp_A(:,:,j),registration2.displacementField,registration2.interpolation);
                end
                links{i}.A=reshape(temp_A,links{i}.imageSize(1)*links{i}.imageSize(2),[]);
                links{i}.A(links{i}.A<10^(-6))=0;
                if ~isempty(links{i}.Cn)
                links{i}.Cn=imtranslate(links{i}.Cn,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                links{i}.Cn=deformation(links{i}.Cn,registration2.displacementField,registration2.interpolation);
                end
                display_progress_bar(i/length(links)*100,false);
            end
        else
            total=0;
            
            for i=1:length(links)
                if isequal(registration_template,'correlation')
                    template_curr=links{i}.Cn;
                    template_curr(template_curr<.9)=0;
                    template_prev=neurons{i}.Cn;
                    template_prev(template_prev<.9)=0;
                else
                    template_curr=max(reshape(links{i}.A./max(links{i}.A,[],1),links{1}.imageSize(1),links{1}.imageSize(2),[]),[],3);
                    template_prev=max(reshape(neurons{i}.A./max(neurons{i}.A,[],1),links{1}.imageSize(1),links{1}.imageSize(2),[]),[],3);
                end
                
                registration1=registration2d(template_prev,template_curr,'transformationModel','translation');
                template=imtranslate(template_curr,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                registration2=registration2d(template_prev,template,'transformationModel',registration_method);
                template_links{i}=deformation(template,registration2.displacementField,registration2.interpolation);
                temp_A=reshape(links{i}.A,links{i}.imageSize(1),links{i}.imageSize(2),[]);
                parfor j=1:size(neurons{i}.A,2)
                    temp_A(:,:,j)=imtranslate(temp_A(:,:,j),-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                    temp_A(:,:,j)=deformation(temp_A(:,:,j),registration2.displacementField,registration2.interpolation);
                end
                links{i}.A=reshape(temp_A,links{i}.imageSize(1)*links{i}.imageSize(2),[]);
                %non-rigid registration can cause spurious non zero values
                links{i}.A(links{i}.A<10^(-6))=0;
                if ~isempty(links{i}.Cn)
                links{i}.Cn=imtranslate(links{i}.Cn,-1*[registration1.transformationMatrix(1,3), registration1.transformationMatrix(2,3)],'FillValues',0);
                links{i}.Cn=deformation(links{i}.Cn,registration2.displacementField,registration2.interpolation);
                end
                total=total+1;
                display_progress_bar(total/length(links)*100,false);
            end
            
        end
        
        display_progress_bar(' Completed',false);
        display_progress_bar('terminating',true)
        
        
    end
    
    disp('Updating Centroids')
    %This will not update if you replace for with parfor. Someone smarter
    %than me could probably tell you why, probably has something to do with
    %class methods in parfor loops
    for i=1:length(neurons)
        neurons{i}.updateCentroid();
    end
    if ~isempty(links)
        for i=1:length(links)
            links{i}.updateCentroid();
        end
    end
end

%% Construct Distance Metrics Between Sessions and Conduct Registration
if ~isempty(links)
    correlation_matrices=construct_correlation_matrix_graph_search(neurons,links,overlap,max_dist);
else
    for i=1:2*length(neurons)
        correlation_matrices{i}=[];
    end
end


%Insert any new distance metrics below, and add them to the
%distance_metrics cell

%Weights order: correlation, centroid_dist,overlap,KL,SNR,decay
%overlap similarity
if weights(3)>0
    overlap_matrices=construct_overlap_matrix_graph_search(neurons);
else
    overlap_matrix=cell(length(neurons)-1,length(neurons));
end

%centroid distance (this metric is required, though  the weight may be 0)
distance_matrices=construct_distance_matrix_graph_search(neurons);





distance_links=cell(1,length(neurons)*2-2);
for i=1:length(neurons)-1
    if ~isempty(links)
        distance_links{2*i-1}=pdist2(neurons{i}.centroid,links{i}.centroid);
        distance_links{2*i}=pdist2(links{i}.centroid,neurons{i+1}.centroid);
    else
        distance_links{2*i-1}=[];
        distance_links{2*i}=[];
    end
end

%JS distance
if weights(4)>0
    JS_matrices=construct_JS_matrix_graph_search(neurons,max_dist);
else
    JS_matrices=cell(length(neurons)-1,length(neurons));
end
if weights(5)>0
for i=1:length(neurons)
    for j=i:length(neurons)
        SNR_dist{i,j}=pdist2(neurons{i}.SNR,neurons{j}.SNR)./((repmat(neurons{i}.SNR,1,size(neurons{j}.SNR,1))+repmat(neurons{j}.SNR',size(neurons{i}.SNR,1),1))/2);
    end
end
else
    SNR_dist=cell(length(neurons),length(neurons));
end
if weights(6)>0
for i=1:length(neurons)
    for j=i:length(neurons)
        decay_dist{i,j}=pdist2(neurons{i}.P.kernel_pars,neurons{j}.P.kernel_pars)./((repmat(neurons{i}.P.kernel_pars,1,size(neurons{j}.P.kernel_pars,1))+repmat(neurons{j}.P.kernel_pars',size(neurons{i}.P.kernel_pars,1),1))/2);
    end
end
else
    decay_dist=cell(length(neurons),length(neurons));
end


%Add additional metrics to the end of this cell. Make sure distance
%matrices is the first element of the cell.
for i=1:size(distance_matrices,1)
    for j=1:size(distance_matrices,2)
        distance_metrics{i,j}={distance_matrices{i,j},overlap_matrices{i,j},JS_matrices{i,j},SNR_dist{i,j},decay_dist{i,j}};
        ind=find(weights(3:end)==0);
        distance_metrics{i,j}(ind+1)=[];
    end
end



% cell of strings 'low', and 'high', indicating whether the distance
% similarity prefers low or high values (centroid distance: low, overlap:
% high) 
% First element should be low, corresponding to centroid distance.
similarity_pref={'low','high','low','low','low'};
similarity_pref(ind+1)=[];

disp('Beginning Cell Tracking')

%vector of weights for total metrics (including correlation) If no links are used,
%set first elements of this vector to 0.
weights=weights/sum(weights);
min_num_neighbors=false;
[Path_Storage,aligned_probabilities]=compute_aligned_pairwise_probability(correlation_matrices,distance_links,distance_metrics,...
    similarity_pref,weights,probability_assignment_method,max_dist,max_gap,min_prob,single_corr,corr_thresh,use_spat,min_num_neighbors,chain_prob,binary_corr);
%neurons=neurons1;
%links=links1;
%% Construct Neuron Throughout Sessions Using Extracted Registration

data_shape=neurons{1}.imageSize;



neuron=Sources2D;

for i=1:size(Path_Storage,1)
    
    A=zeros(size(neurons{1}.A(:,1)));
    C={};
    S={};
    C_raw={};
    for k=1:size(Path_Storage,2)
        if ~isnan(Path_Storage(i,k))
            C{k}=neurons{k}.C(Path_Storage(i,k),:);
            C_raw{k}=neurons{k}.C_raw(Path_Storage(i,k),:);
            S{k}=neurons{k}.S(Path_Storage(i,k),:);
            A=A+neurons{k}.A(:,Path_Storage(i,k));
        else
            C{k}=zeros(1,size(neurons{k}.C,2));
            C_raw{k}=zeros(1,size(neurons{k}.C,2));
            S{k}=zeros(1,size(neurons{k}.C,2));
        end
    end
    
    
    A=A/size(Path_Storage,2);
    neuron.A=horzcat(neuron.A,A);
    neuron.C=vertcat(neuron.C,horzcat(C{:}));
    neuron.C_raw=vertcat(neuron.C_raw,horzcat(C_raw{:}));
    neuron.S=vertcat(neuron.S,horzcat(S{:}));
    centroid=calculateCentroid(A,data_shape(1),data_shape(2));
    neuron.centroid=vertcat(neuron.centroid,centroid);
    decays=[];
    for q=1:size(Path_Storage,2)
        decays=[decays,neurons{q}.P.kernel_pars(Path_Storage(i,q))];
    end
    neuron.P.kernel_pars(i,1)=mean(decays);
  
end
neuron.Cn=neurons{base}.Cn;
neuron.imageSize=neurons{base}.imageSize;
neuron.updateCentroid();
neuron=calc_snr(neuron);
neuron.probabilities=aligned_probabilities;

neuron.existant_indices=Path_Storage;
try
    neuron.options=neurons{1}.options;
end
%Eliminate remaining neurons falling below chain_prob
neuron.delete(neuron.probabilities<chain_prob);
