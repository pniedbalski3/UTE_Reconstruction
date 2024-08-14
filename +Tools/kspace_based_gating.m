function binning = kspace_based_gating(fid,coils,NBins)
if mod(NBins,2) ~=0
    disp(['Number of Bins must be even, using NBins + 1 = ' num2str(NBins+1) ' bins.']);
    NBins = NBins+1;
end

%We want this pretty smooth - TR in the ballbark of 5 ms. Resp in the
%ballpark of 4 s/cycle. smooth with a window of 200 points
gate_k0 = smooth(sqrt(sum(abs(fid(1,:,coils)).^2,3)),200);
gate_k0(1:100) = [];
gate_k0((end-99):end) = [];

for i = 1:(length(gate_k0)/5000+1)
    if i*5000 < length(gate_k0)
        gate_k0(((i-1)*5000+1):(i*5000)) = gate_k0(((i-1)*5000+1):(i*5000))/mean(gate_k0(((i-1)*5000+1):(i*5000)));
    else
        gate_k0(((i-1)*5000+1):end) = gate_k0(((i-1)*5000+1):end)/mean(gate_k0(((i-1)*5000+1):end));
    end
end
% Assume we'll have top and bottom bins + 2 bins for each intermediate
% respiratory phase
nclusters = (NBins-2)/2 + 2;

% use k-means clustering for a really easy way to get a coarse binning
binning = kmeans(gate_k0,nclusters);

%Need to sort clusters for easier handling
clust_means = zeros(nclusters,1);
for i = 1:nclusters
    clust_means(i) = mean(gate_k0(binning == i));
end
[~,sortind] = sort(clust_means);

for i = 1:length(binning)
    newbin = find(sortind == binning(i));
    binning(i) = newbin;
end


%% Now, need to divide each of the intermediate bins to account for during inspiration vs. during expiration.

dgate_k0 = smooth(smooth(diff(gate_k0),200),200);

dgate_k0(end+1) = dgate_k0(end);

for i = 1:length(gate_k0)
    if binning(i) == 1 || binning(i) == nclusters
        continue
    else
        if dgate_k0(i) < 0
            binning(i) = NBins-binning(i)+2;
        end
    end
end

% start_ind = 1;
% mybin = binning(1);
% for i = 1:length(gate_k0)
%     %If we're at the highest or lowest bin, don't separate into two.
%     if mybin == 1 || mybin == nclusters
%         mybin = binning(i);
%         start_ind = i;
%         continue;
%     end
%     if binning(i) ~= mybin
%         stop_ind = i-1;
%         diff_val = mean(dgate_k0(start_ind:stop_ind));
%         if diff_val < 0 
%             binning(start_ind:stop_ind) = NBins-mybin+2;
%         end
%         start_ind = i;
%         mybin = binning(i);
%     end
% end

% Remove points from inspiration that are too far off - Simple way is to
% look at Expiration (which should be pretty stable) and make sure
% inspiration and expiration have the same bin width. I'll have to check
% some other people to see if this works in all cases
exp_range = max(gate_k0(binning == NBins/2+1)) - min(gate_k0(binning == NBins/2+1));

binning(gate_k0 < max(gate_k0(binning ==1) - exp_range)) = NaN;

%% Plot the final binning.
figure('Name','k0-based_Binning');
set(gcf,'Position',[77 452 1660 420]);
hold on
xax = 1:length(binning);
for i = 1:NBins
    plot(xax(binning==i),gate_k0(binning==i),'*');
end




