function moving_proj_temp=autoAlignment_adj(vid_paths,fixed_index,projection_paths)
%Modifies vid_files listed in vid_paths to match the first file in vid_path

%For example, use this code to consecutively align videos in a folder, or align all videos to a baseline.

%n -no, y - yes, m - manual intervention


close all
if ~exist('fixed_index','var')||isempty(fixed_index)
    fixed_index=1;
end
for i=1:length(vid_paths)
    Y=load(vid_paths{i});
    try
    vids{i}=Y.Y;
    catch
        vids{i}=Y.current;
    end
    
    
end
try
for i=1:length(projection_paths)
    project=struct2cell(load(projection_paths{i}));
    projections{i}=project{1};
end
end
if ~exist('projections','var')||isempty(projections)
    projections={};
end

aligned={};
aligned{fixed_index}=vids{fixed_index};
nonempty=1;
while  nonempty<length(vid_paths)
    while true
        [~,indices]=number_nonempty(aligned);
        indices_string="";
        for i=1:length(indices)
            indices_string=strcat(indices_string,' ',string(indices(i)));
        end
        base_selection=input(strcat('select video index to be used as base, ',indices_string,': '),'s');
        try
            base_selection=str2num(base_selection);
            if any(indices==base_selection)
                break
            end
        end
    end
    fixed=double(aligned{base_selection});
    unaligned=setdiff(1:length(vid_paths),indices);
    while true
       
        indices_string="";
        for i=1:length(unaligned)
            indices_string=strcat(indices_string,' ',string(unaligned(i)));
        end
        alter_selection=input(strcat('select video index to alter, ',indices_string,': '),'s');
        try
            alter_selection=str2num(alter_selection);
            if any(unaligned==alter_selection)
                break
            end
        end
    end
    moving=double(vids{alter_selection});
    
    
    
    
    
    
    
    
%     
%     
%     light=[];
%     fix_back=zeros(size(fixed));
%     parfor k=1:size(fixed,3)
%         fix_back(:,:,k)=medfilt2(remove_template_background(fixed(:,:,k)));
%         light(k)=sum(fixed(:,:,k),'all');
%     end
%     prc=prctile(light,60);
%     indices=find(light<prc);
%     
%     fixed_proj1=std(fix_back(:,:,indices),[],3);
%     fixed_proj2=max(fix_back(:,:,indices),[],3);
%     fixed_proj=(fixed_proj1+fixed_proj2)/2;
%     light=[];
%     mov_back=zeros(size(moving));
%     parfor k=1:size(moving,3)
%         mov_back(:,:,k)=medfilt2(remove_template_background(moving(:,:,k)));
%         light(k)=sum(moving(:,:,k),'all');
%     end
%     prc=prctile(light,60);
%     indices=find(light<prc);
%     
%     
%     moving_proj1=std(mov_back(:,:,indices),[],3);
%     moving_proj2=max(mov_back(:,:,indices),[],3);
%     moving_proj=(moving_proj1+moving_proj2)/2;
    if length(projections)<base_selection||isempty(projections{base_selection})
        projections{base_selection}=construct_template(fixed);
    end
    
    if length(projections)<alter_selection||isempty(projections{alter_selection})
        projections{alter_selection}=construct_template(moving);
    end
    fixed_proj=projections{base_selection};
    moving_proj=projections{alter_selection};
    fixed_proj(fixed_proj<.9)=0;
    moving_proj(moving_proj<.9)=0;
    registration=registration2d(fixed_proj,moving_proj);
    registration=registration2d(fixed_proj,moving_proj);
    
    fixed_proj(fixed_proj<.9)=0;
    moving_proj(moving_proj<.9)=0;
    moving_proj_temp=deformation(moving_proj,registration.displacementField...
        ,registration.interpolation);
    
    
    figure('name','moving','Position', [10 10 200 200])
    imagesc(moving_proj_temp)
    colormap gray
    daspect([1,1,1])
    figure('name','fixed','Position', [10 10 200 200])
    imagesc(fixed_proj)
    colormap gray
    daspect([1,1,1])
    
    %m denotes manual feature selection
    keep=input('Is registration acceptable? y/n/m ','s');
    %keep='y';
    while isequal(keep,'m')
        
        fixed_proj(fixed_proj<.9)=0;
        moving_proj(moving_proj<.9)=0;
        [fixed_mask,moving_mask]=manual_ROI_selection(fixed_proj,moving_proj);
        for i=1:length(fixed_mask)
            s=regionprops(fixed_mask{i},'centroid');
            centroid_1(i,:)=s.Centroid;
            s=regionprops(moving_mask{i},'centroid');
            centroid_2(i,:)=s.Centroid;
        end
        diff=centroid_2-centroid_1;
        diff=mean(diff,1);
        
        moving_proj_temp=imtranslate(moving_proj,-1*diff);
        
        
        registration=registration2d(fixed_proj,moving_proj_temp,'transformationModel','translation');
        moving_proj_temp=deformation(moving_proj_temp,registration.displacementField...
            ,registration.interpolation);
        figure('name','moving','Position', [10 10 200 200])
        imagesc(moving_proj_temp)
        colormap gray
        daspect([1,1,1])
        figure('name','fixed','Position', [10 10 200 200])
        imagesc(fixed_proj)
        colormap gray
        daspect([1,1,1])
        keep=input('Is registration acceptable? y/n/m ','s');
    end
    
    
    
    
    
    if isequal(lower(keep),'y')
        if exist('diff','var');
            parfor k=1:size(moving,3);
                moving(:,:,k)=imtranslate(moving(:,:,k),-1*diff);
            end
            clear diff
        end
        
        parfor k=1:size(moving,3)
            moving(:,:,k) = uint8(deformation(double(moving(:,:,k)),...
                registration.displacementField,registration.interpolation));
            
        end
        projections{alter_selection}=moving_proj_temp;
  
        disp('registration accepted')
        aligned{alter_selection}=moving;
  
    end
    
    saver=input('save current aligned files, and end session? y/n ','s');
    if isequal(saver,'y')
        saved_files={};
        for i=1:length(aligned)
            Y=uint8(aligned{i});
            if ~isempty(Y)
                Ysiz=size(Y);
                save(vid_paths{i},'Y','Ysiz','-v7.3');
                try
                template=projections{i};
                save(projection_paths{i},'template')
                end
                saved_files{i}=vid_paths{i};
            end
        end
        disp('saved files')
        disp(saved_files)
        return
    end
    
    [nonempty,indices]=number_nonempty(aligned);
    close all
end
for i=1:length(projections)
    figure()
    imagesc(projections{i});
    colormap gray
    daspect([1,1,1])
end



saved_files={};
for i=1:length(aligned)
    Y=uint8(aligned{i});
    Ysiz=size(Y);
    save(vid_paths{i},'Y','Ysiz','-v7.3');
    saved_files{i}=vid_paths{i};
end
try
for i=1:length(aligned)
    template=projections{i};
    save(projection_paths{i},'template')
end
end



disp('saved files')
disp(saved_files)




        
        
        