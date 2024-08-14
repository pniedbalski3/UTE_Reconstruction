function Im = simple_bart_recon(fid,traj,ImSize,FOV,weights,tmp_path)
%%
    % Generate density compensation weights
    nIter = 10;
    myDCF = Recon.get_DCF(traj,ImSize,nIter);

    barttraj = traj*FOV*0.8375; %Scale for discrepancy between the Image size listed in the raw data and the actual image size acquired
    bartdcf = reshape(myDCF,[1 size(myDCF)]);
    bartfid = reshape(fid,[1 size(fid)]);
    
    bartfid = bartfid(:,:,1:length(weights),:);
    barttraj = barttraj(:,:,1:length(weights));
    bartdcf = bartdcf(:,:,1:length(weights));

    tmpstr = [num2str(ImSize) ':' num2str(ImSize) ':' num2str(ImSize) ' '];
    % Scale data by dcf
    datac = bart('fmac',bartfid,bartdcf);
    %Base Recon
    imgL = bart(['nufft -a -d ' tmpstr],barttraj,datac);
    % FFT
    ksp = bart('fft 7',imgL);
    % Coil Combination - Use 8 coils to make sure we have adequate memory
    ccMatrix = bart('cc -M',ksp);
    data1 = bart('ccapply -p 8',bartfid,ccMatrix);
    ksp1 = bart('ccapply -p 8',ksp,ccMatrix);
    % Estimate coil sensitivies from k-spac>getNextLine (line 38)
    mapsL = bart('caldir 24',ksp1);

    filestr = fullfile(tmp_path,['tmpdcfdata' num2str(round(rand(1)*1000))]);

    %Soft gating step
    sg_dcf = bartdcf;
    if length(weights) >1
        for i = 1:length(weights)
            sg_dcf(:,:,i) = sqrt(bartdcf(:,:,i))*weights(i);
        end
    else
        sg_dcf = sqrt(bartdcf);
    end

    writecfl(filestr,sg_dcf);
    
    filestr = strrep(filestr,'\','/');
    filestrWSL = wslPathCorrection(filestr);
    %Put in quotes just in case windows path string has spaces
    filestrWSL2 = ['"' filestrWSL '"'];
    %Final Recon
    Im = bart(['pics -C 50 -i 80 -R T:7:0:0.001 -p ' filestrWSL2 ' -t'],barttraj,data1,mapsL);