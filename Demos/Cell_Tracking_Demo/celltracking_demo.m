%Parameters are detailed in cellTracking_SCOUT
clear all; clc
load('neurons')
load('links')
neuron=cellTracking_SCOUT(neurons(1:2),'links',links(1:2),'overlap',1950,'register_sessions',true,'weights',[5,5,5,5,0,0],'probability_assignment_method','Kmeans','registration_template','spatial','max_dist',14,'registration_method','non-rigid','chain_prob',.9,'corr_thresh',0,'min_prob',.5);