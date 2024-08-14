function [slice,dim,coils] = ONE_select_slice_coils(fid,traj,Orig_ImSize,Desired_Size,N_short)

%In this function, just recon a small portion of the data for every coil
%element so that we can select which slice best shows diaphragm and which
%coil elements are nearest that region.

if nargin<5
    N_short = 200;
end
%%
mytraj = traj(:,:,1:N_short);
myfid = fid(:,1:N_short,:);

[myfid,mytraj] = Tools.change_mat_size(myfid,mytraj,Desired_Size,Orig_ImSize);

%% Reconstruct
%Make Image to hold everything
Image = zeros(Desired_Size,Desired_Size,Desired_Size,size(myfid,3)+1);

nIter = 10;
myDCF = Recon.get_DCF(mytraj,Desired_Size,nIter);

parfor i = 1:size(myfid,3)
    Image(:,:,:,i+1) = Recon.mem_eff_recon(Desired_Size,myfid(:,:,i),mytraj,myDCF,i,size(myfid,3));
end
Image(:,:,:,1) = Tools.soscoilcombine(Image); 

cent = size(Image,1)/2;

figure('Name','Pick Slice')
set(gcf,"Position",[959 42 863 954]);
a = 6;b = 6;
tiledlayout(a,b,"TileSpacing",'none')
dim = 2;
slice = 1;
while slice == 1
    moveind = -(floor(a*b/2)):(floor(a*b/2));
    for i = 1:(a*b)
        if dim == 2
            myslice = Image(:,cent+moveind(i),:,1);
        elseif dim == 1
            myslice = Image(cent+moveind(i),:,:,1);
        elseif dim == 3
            myslice = Image(:,:,cent+moveind(i),1);
        end
        nexttile;
        imagesc(squeeze(abs(myslice)/max(abs(myslice(:)))));
        colormap(gray);axis square;axis off;clim([0 0.6]);
        title(num2str(cent+moveind(i)))
    end
    prompt = {'What slice best shows diaphragm?','If Dimension not correct, try 1,2,or 3'};
    dlgtitle = 'Diaphragm Slice Selection';
    fieldsize = [1 40];
    definput = {'1','2'};
    opts.Interpreter = 'tex';
    answer = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    slice = str2double(answer{1});
    dim = str2double(answer{2});
    close;
end

figure('Name','Pick Coils')
set(gcf,'Position',[583 42 1230 954])
for i = 1:size(myfid,3)
    if dim == 2
        myslice = Image(:,slice,:,i+1);
    elseif dim == 1
        myslice = Image(slice,:,:,i+1);
    elseif dim == 3
        myslice = Image(:,:,slice,i+1);
    end
    nexttile;
    imagesc(squeeze(abs(myslice)/max(abs(myslice(:)))));
    colormap(gray);axis square;axis off;clim([0 0.6]);
    title(num2str(i))
end
prompt = {'What coils best pick up the diaphragm (pick up to 2. Separate by commas if necessary)?'};
dlgtitle = 'Coil Selection';
fieldsize = [1 40];
definput = {'7,8'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
coils = str2num(answer{1});
close;

