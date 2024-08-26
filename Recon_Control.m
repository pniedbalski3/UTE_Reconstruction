%--------------------------------------------------------------------------
% This script walks through the steps required to reconstruct UTE MRI data acquired during free breathing.
% Options are provide for k0-based gating and image-based gating as well as
% hard and soft gating. 
% A protocol describing the use of this code for the generation of UTE MRI
% has been submitted to the Journal of Visualized Experiments (JoVE)
%
%
% Author: Peter J. Niedbalski
% pniedbalski@kumc.edu
%--------------------------------------------------------------------------

%% Clean Workspace
clear;clc;close all;

%% Step 1: Load UTE Data

%function to load radial imaging data for an in-house 3D UTE MRI sequence
%built for a Siemens 3T MRI system

[fid,traj,ImSize,FOV,mypath] = DataImport.load_radial_data();

% The rest of the code expects:
% fid: complex array of size Npts x Nprojections x nCoilElements
% traj: array of size 3 x Npts x Nprojections
% ImSize: Scalar
% FOV: Scalar
% mypath: string - location of the data and where images will be written
% out by default

Orig_ImSize = 320;

%Trajectory delay correction in units of dwell time (i.e. shift acquisition
%by 4 points to account for slow turn on of gradients. This may need to be
%adjusted for a given sequence or scanner
traj = DataImport.traj_delay_correction(traj,4);

%% Step 2: Identify coil elements closest to diaphragm and slice/dimension for image-based gating
%Pass fids, trajectories, the collected image size (will be subsampled to
%96x96x96 matrix), and the number of projections to use (200 seems to be
%adequate - makes for a fast recon)
Desired_Size = 96;
NProj = 200;

[slice,dim,coils] = ReconSteps.ONE_select_slice_coils(fid,traj,Orig_ImSize,Desired_Size,NProj);
Sliding_Window = NProj;

%% Step 3: Sliding Window Reconstruction for Image-based gating
All_Im = ReconSteps.TWO_sliding_window_recon(fid,traj,Orig_ImSize,slice,dim,coils,Desired_Size,Sliding_Window);
All_Im = abs(All_Im);

%% Step 4: Select line over diaphragm and use 1-D navigator to bin respiratory states
good = 0;
while good == 0
    figure('Name','Select indices of line over diaphram and press enter')
    set(gcf,'Position',[1020 223 834 688]);
    imagesc(squeeze(abs(All_Im(:,:,1))))
    colormap(gray)
    axis square;
    axis off;
    [xpoints,ypoints] = ginput;
    close;
    xpoints = round(xpoints);
    ypoints = round(ypoints);
    if length(xpoints) == 2
        diffx = xpoints(2) - xpoints(1);
        diffy = ypoints(2) - ypoints(1);
        if diffx > diffy
            ypoints(1) = round(mean(ypoints));
            ypoints(2) = ypoints(1);
        elseif diffy > diffx
            xpoints(1) = round(mean(xpoints));
            xpoints(2) = xpoints(1);
        end
        good = 1;
    end
end

diaphragm_pos = ReconSteps.THREE_diaphragm_motion(All_Im,xpoints,ypoints);

%% Convert Respiratory motion to actual gating indices
%From the diaphragm positions found in the previous step, create an array
%holding projection indices for each diaphram position.
my_index = ReconSteps.FOUR_get_gated_indices(diaphragm_pos,Sliding_Window,fid);

%% Generate soft-gating weights
%Assume end expiration is the state with the most binned points.
numbinned = sum(my_index,2);
[~,final_gate] = max(numbinned);

%Simple soft-gating
myweights = Tools.image_based_soft_gating(my_index,final_gate);
%For hard-gating instead of soft gating:
%myweights = 1;

%Sanity check to look at the gating indices.
figure('Name','Check_Gating_Indices')
imagesc(my_index);
set(gcf,'Position',[37 475 1797 420]);
drawnow

%% Reconstruction
% Data can't be passed to WSL from a secure network location, so pass to a
% temp folder

tmpdir = 'C:\Users\pniedbalski\OneDrive - University of Kansas Medical Center\Documents\tmp_bart_analysis';
cd(tmpdir)

Image = Recon.simple_bart_recon(fid,traj,Orig_ImSize,FOV,myweights,tmpdir);

%% Write data to file (default to same location as raw data)
cd(mypath)
niftiwrite(abs(Image),'ImBased_SoftGate_EndExp_Recon_Image','Compressed',true);

%% Alternatively, perform k0-based gating
%Very important that we're at steady state magnetization: Delete a couple
%thousand points:
tmpfid = fid;
tmpfid(:,1:3000,:) = [];
%In the binning function, we'll be deleting 200 points from the beginning
%and end (to account for issues with smoothing), so get rid of those in the
%trajectory array
tmptraj = traj;
tmptraj(:,:,[1:3100 (end-99):end]) = [];

%Number of respiratory bins to use for k-space gating
NBins = 6;
binning = Tools.kspace_based_gating(tmpfid,coils,NBins);
tmpfid(:,[1:100 (end-99):end],:) = [];

%% Generate soft-gating weights
tmpfid(:,isnan(binning),:) = [];
tmptraj(:,:,isnan(binning)) = [];

kweights = Tools.kspace_based_soft_gating(binning,4);
kweights(isnan(binning)) = [];
%for hard gating instead of soft-gating:
%kweights = 1;

%% Perform Image reconstruction
cd(tmpdir)
Image = Recon.simple_bart_recon(tmpfid,tmptraj,Orig_ImSize,FOV,kweights,tmpdir);

%% Write data to file (default to same location as raw data)
cd(mypath)
niftiwrite(abs(Image),'k0Based_SoftGate_EndExp_Recon_Image','Compressed',true);