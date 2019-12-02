function full_neuron=full_demo_endoscope(filename,indices,JS,min_pnr,JS_constraint,gSiz)
%vid_file: (m x n x num_frames) matrix containing video data

% indices: specify frame indices you wish to perform extraction on, leave as [] for all frames

% JS constraint: percentile, between 0 and 1. Anything outside the
% percentile range will be eliminated
% min_pnr: min allowable peak-to-noise ratio. smaller values pick up more neurons, but cause more false positives. typical values, (5:15)





%% clear workspace

global  d1 d2 numFrame sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights; %#ok<NUSED> % global variables, don't change them manually



ssub=2;
tsub=2;

if ~exist('JS','var')||isempty('JS')
   JS=0;
end

if ~exist('JS_constraint','var')||isempty(JS_constraint)
    JS_constraint='bound'; %Options 'bound' (JS variable is the actual bound on the maximum JS value) 
                      %prc (JS variable indicates a percentile range.
                      %Anything outside the range is deleted. Typically .95
                      %or .9
end
divergence_type='JS'; %JS or JS, type of divergence used to calculate differences between probability distributions.

%% select data and map it to the RAM
nam = filename;

cnmfe_choose_data;

%% create Source2D class object for storing results and parameters
Fs = 15;             % frame rate
        % spatial downsampling factor
         % temporal downsampling factor
gSig = 3;           % width of the gaussian kernel, which can approximates the average neuron shape
if ~exist('gSiz','var')||isempty(gSiz)
    gSiz = 24;          % maximum diameter of neurons in the image plane. larger values are preferred.
    
end
gSiz=gSiz/2;
gSizMin=7;          % minimum diameter of neurons in the image plane.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron_full.Fs = Fs;         % frame rate

% with dendrites or not 
with_dendrites = false;
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 20;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 3;
end
%neuron.JS=JS;

%neuron_full.JS=JS;
merge_thr = [0.5, 0.75, -1];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
dmin =[1,15];
%dmin=1;
%% options for running deconvolution 
neuron_full.options.deconv_flag = true; 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

%% downsample data for fast and better initialization
if ~isempty(indices)
    sframe=indices(1);						% user input: first frame to read (optional, default:1)
    num2read= indices(2)-indices(1);             % user input: how many frames to read   (optional, default: until the end)

else
    sframe=1;
    num2read=Ysiz(3);
end
tic;
cnmfe_load_data;
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);



Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix

%% compute correlation image and peak-to-noise ratio image.
%cnmfe_show_corr_pnr;    % this step is not necessary, but it can give you some...
                        % hints on parameter selection, e.g., min_corr & min_pnr

%% initialization of A, C
% parameters
debug_on = false;   % visualize the initialization procedue. 
save_avi = false;   %save the initialization procedure as an avi movie. 
patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large 
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = 0.8;     % minimum local correlation for a seeding pixel
if ~exist('min_pnr','var')||isempty('min_pnr')
    min_pnr = 8;       % minimum peak-to-noise ratio for a seeding pixel
end
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron_full.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron.options.nk = 5;  % number of knots for detrending 
%neuron.imageSize=[ceil(d1/ssub),ceil(d2/ssub)];
% greedy method for initialization
tic;

[center,Cn,pnr]=neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
           
         

fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

% show results
figure;
imagesc(Cn, [0.1, 0.95]);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;

% sort neurons
neuron.orderROIs('snr');
neuron_init = neuron.copy();

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = ssub;    % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 5;       % this variable has different meanings: 
                    %1) udpate_spatial_method=='hals' or 'hals_thresh',
                    %then Nspatial is the maximum iteration 
                    %2) update_spatial_method== 'nnls', it is the maximum
                    %number of neurons overlapping at one pixel 
               
% parameters for running iteratiosn 
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 1;        % maximum number of iterations 
miter = 1; 
while miter <= maxIter
        % if JS>0
           
        %    [~,JS_score]=spatial_filter(neuron,JS,[ceil(d1/ssub),ceil(d2/ssub)],JS_constraint,true,gSiz,gSizMin,true,divergence_type);
         
        %end

    %% merge neurons
    neuron.remove_false_positives();
    try
    cnmfe_quick_merge;              % run neuron merges
      cnmfe_merge_neighbors;          % merge neurons if two neurons' peak pixels are too close 
    end
    %% udpate background 
    % estimate the background
    tic;
    cnmfe_update_BG;
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
       
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial & temporal components
    tic;
    for zz=1:2
         neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();    % run this line if neuron shapes are circular 
  
        %neuron=split_neurons(neuron,[ceil(d1/ssub),ceil(d2/ssub)],dmin,gSizMin);
        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
         neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();    % run this line if neuron shapes are circular 
  
  
        if JS>0
          
           [~,JS_score]=spatial_filter(neuron,JS,[ceil(d1/ssub),ceil(d2/ssub)],JS_constraint,true,gSiz,gSizMin,true,'threshold',divergence_type);
         
        end
        
       
        %temporal
        neuron.updateTemporal_endoscope(Ysignal);
        try
        cnmfe_quick_merge;              % run neuron merges
        end
        
        neuron.remove_false_positives();
        try
        cnmfe_merge_neighbors; 
        end
%         
       
    end
     
     %neuron=filter_few_spikes(neuron,20,.025);
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
       K=200;
     %% pick neurons from the residual (cell 4).
     for resiter=1

         
         
          
      if miter==1
         seed_method = 'auto'; 
         [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, seed_method, debug_on,K); % method can be either 'auto' or 'manual'
     end
     
     %if JS>0
           
     %       [~,JS_score]=spatial_filter(neuron,JS,[ceil(d1/ssub),ceil(d2/ssub)],JS_constraint,true,gSiz,gSizMin,false,divergence_type);
         
     %end
        
     try
     cnmfe_quick_merge;              % run neuron merges
    cnmfe_merge_neighbors;          % merge neurons if two neurons' peak pixels are too close 
     end
    %% udpate background 
    % estimate the background
    tic;
    
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial & temporal components
    tic;
    for zz=1:2
        
        
        cnmfe_update_BG;
        fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
        
       
        %temporal
        try
        neuron.updateTemporal_endoscope(Ysignal);
        end
   
        try
        cnmfe_quick_merge;              % run neuron merges
        end
        %spatial
     
  
         neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();    % run this line if neuron shapes are circular 
  
        %neuron=split_neurons(neuron,[ceil(d1/ssub),ceil(d2/ssub)],dmin,gSizMin);

        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
         neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();    % run this line if neuron shapes are circular 
  
  
        if JS>0
          
           [~,JS_score]=spatial_filter(neuron,JS,[ceil(d1/ssub),ceil(d2/ssub)],JS_constraint,true,gSiz,gSizMin,true,'threshold',divergence_type);
         
        end
 
%         % post process the spatial components (you can run either of these two operations, or both of them)
             neuron.remove_false_positives();
        try
        cnmfe_merge_neighbors; 
        end
        % stop the iteration when neuron numbers are unchanged. 
        %if isempty(merged_ROI)
        %    break;
        %end
    end
     
    neuron.remove_false_positives(); 
        end
      neuron.compactSpatial(); 
       neuron.compactSpatial(); 
       miter=miter+1;
end
neuron_downsample=neuron.copy();     
if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    neuron.imageSize=[d1,d2];
    cnmfe_full;
    neuron_full = neuron.copy();
end    
    
%% sort neurons
[~, srt] = sort(max(neuron.C, [], 2), 'descend'); % C: row is the number of neuron, column is the number of frames
neuron.orderROIs(srt);
neuron_init = neuron.copy();
%  if JS>0
%            
%            [~,JS_score]=spatial_filter(neuron,JS,[ceil(d1),ceil(d2)],JS_constraint,false,gSiz*2,gSizMin*2,false,divergence_type);
%          
%          
%        end
      

gSiz=gSiz*2;
%spatial
neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
 
%temporal
try
neuron.updateTemporal_endoscope(Ysignal);
end
 
%% delete some neurons and run CNMF-E iteration 
neuron.orderROIs('decay_time');  % you can also use {'snr', 'mean', 'decay_time'} 

tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
%update spatial & temporal components
tic;
%neuron=smooth_spatial(neuron);
%neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
for zz=1:2
        %if zz==1
        neuron.updateSpatial_endoscope(Ysignal,Nspatial,update_spatial_method);
       % neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
       %  neuron.compactSpatial(); 
        
        %neuron=split_neurons(neuron,[d1,d2],dmin,gSizMin);
    
        %neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
         %neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        %neuron.compactSpatial();    % run this line if neuron shapes are circular 
  
  
        if JS>0
          
           [~,JS_score]=spatial_filter(neuron,JS,[d1,d2],JS_constraint,true,gSiz,gSizMin,true,'threshold',divergence_type);
         
        end
 
        
    %temporal
   try
    neuron.updateTemporal_endoscope(Ysignal);
    
   end
    try
    cnmfe_quick_merge;              % run neuron merges
    end
    %spatial
   
       neuron.remove_false_positives();
    try
        cnmfe_merge_neighbors; 
    end
        %neuron=filter_few_spikes(neuron,20,.025);

end
 
neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
neuron.compactSpatial();    % run this line if neuron shapes are circular 

 
if JS>0
          
    
           [~,JS_score]=spatial_filter(neuron,JS,[d1,d2],JS_constraint,true,gSiz,gSizMin,true,'threshold',divergence_type);
         
end
        
%quickMerge(neuron,[.75,.7,-1],gSiz);
%quickMerge(neuron,[.5,.9,-1],gSiz);


neuron.compactSpatial();

neuron.remove_false_positives();

fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

b0 = mean(Y,2)-neuron.A*mean(neuron.C,2); 
Ybg = bsxfun(@plus, Ybg, b0-mean(Ybg, 2)); 
neuron.orderROIs('snr'); 
%

figure;
Cn = imresize(Cn, [d1, d2]);
[neuron.Coor,json_file,centroid] = plot_contours(neuron.A, Cn, 0.8, 1, [], [], 2);
% plot_contours(neuron.A, Cn, 0.8, 0, [], [], 2); % display no numbers
colormap winter;
title('contours of estimated neurons');
%%savejson('jmesh',json_file,'contoursJason8thGsig8Ds3manual2Auto'); %%Suoqin added this to save the image

figure
imagesc(Cn,[min(Cn(:)) max(Cn(:))])
% hold on
% scatter(centroid(:,2),centroid(:,1),'r')
axis tight; axis equal;axis off
colormap winter;

neuron.file = neuron;
neuron.numFrame = numFrame;
neuron.num2read = num2read;
neuron.centroid = centroid;

neuron.Cn = Cn;
[C_df,C_raw_df, Df] = extract_DF_F_endoscope(neuron,Ybg);
neuron.C_df = C_df;neuron.Df = Df;
neuron=calc_snr(neuron);
%indices=find(neuron.SNR<1);
%neuron.delete(indices);
%neuron=filter_SNR(neuron,.1);
neuron.trace = neuron.C.*(max(neuron.A)'*ones(1,size(neuron.C,2)));

close all
%neuron=threshold_C(neuron);
neuron.imageSize=[d1,d2];
neuron.A=full(neuron.A);
%[~,avi_name,~]=fileparts(filename);
%neuron.runMovie(Ysignal, [0, size(Ysignal,3)], true, avi_name);
full_neuron=neuron.copy();

end

