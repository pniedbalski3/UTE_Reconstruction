function All_Im = TWO_sliding_window_recon(fid,traj,Orig_ImSize,slice,dim,coils,Small_size,Sliding_window_size)

%Function to do a sliding window and reconstruct tons! of images. These
%images are 3D, and to save >1000 would be insane. Just save the slice that
%is passed by slice. To save time, only reconstruct the coils closest to
%the diaphragm, which is specified by coils.
if nargin < 8
    Sliding_window_size = 200; %200 seems to work reasonably well.
end

fid = fid(:,:,coils);

if Orig_ImSize ~= Small_size
    [myfid,mytraj] = Tools.change_mat_size(fid,traj,Small_size,Orig_ImSize);
else
    myfid = fid;
    mytraj = traj;
end

%Overlap by 50% - This may be too slow... We'll see.

%get total number of Images in a lazy way
Ind_Ims_Sw = 1:(Sliding_window_size/2):(size(fid,2)-Sliding_window_size);
Ims = zeros(Small_size,Small_size,length(Ind_Ims_Sw));
%tmpIms = zeros(Small_size,Small_size,Small_size,length(coils));
nIter = 10;
parfor i = 1:length(Ind_Ims_Sw)
    disp(['Starting Recon ' num2str(i) ' of ' num2str(length(Ind_Ims_Sw))]);
    tmpfid = myfid(:,Ind_Ims_Sw(i):(Ind_Ims_Sw(i)+Sliding_window_size),:);
    tmptraj = mytraj(:,:,Ind_Ims_Sw(i):(Ind_Ims_Sw(i)+Sliding_window_size));
    tmpDCF = Recon.get_DCF(tmptraj,Small_size,nIter);
    tmpImsA = Recon.mem_eff_recon(Small_size,tmpfid(:,:,1),tmptraj,tmpDCF,1,size(myfid,3),[0 0 0],false);
    if length(coils)>1
        tmpImsA(:,:,:,2) = Recon.mem_eff_recon(Small_size,tmpfid(:,:,2),tmptraj,tmpDCF,2,size(myfid,3),[0 0 0],false);
        tmpIm = Tools.soscoilcombine(tmpImsA);
    else
        tmpIm = tmpImsA;
    end
    if dim==1
        Ims(:,:,i) = squeeze(tmpIm(slice,:,:));
    elseif dim == 2
        Ims(:,:,i) = squeeze(tmpIm(:,slice,:));
    elseif dim == 3
        Ims(:,:,i) = squeeze(tmpIm(:,:,slice));
    else
        error('dim needs to be 1,2, or 3');
    end
    disp(['Done with Recon ' num2str(i) ' of ' num2str(length(Ind_Ims_Sw))]);
end

All_Im = Ims;