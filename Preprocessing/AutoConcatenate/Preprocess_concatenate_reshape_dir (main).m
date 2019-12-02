function []=Preprocess_concatenate_reshape_dir(base_dir,save_dir)
%Takes in the base dir containing all desired files, and recursively
%searches through directories, locating any directory containing files with
%name msCam(i), where i is an integer. The msCam files from each directory
%are then concatenated and reshaped, and the resulting files are saved in
%the save_dir.

%If connection is lost, some files may not be concatenated, the base
%folders will be displayed if this occurs.

%On my linux, subdir required me to navigate to the subdir, and use './' (.\ in windows) as
%the base_dir input


if isequal(save_dir(end),'/')||isequal(save_dir(end),'\')
    save_dir=save_dir(1:end-1);
end



files=subdir(base_dir);
files={files.name};
folders={};
for i=1:length(files)
    [folders{end+1},~,~]=fileparts(files{i});
end
folders=unique(folders);
% for i=length(folders):-1:1
%     [~,name,~]=fileparts(folders{i});
%     if isempty(regexp(name,'m2+'))
%         folders(i)=[];
%     end
% end


parfor i=1:length(folders)
   if ~isempty(strfind(folders{i},'m2')) 
   %Preprocess_concatenate_reshape(folders{i},save_dir);
   concatenate_behav(folders{i},save_dir);
   end
end
