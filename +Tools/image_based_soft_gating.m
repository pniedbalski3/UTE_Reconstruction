function weights = image_based_soft_gating(GateArray,use_ind)

weights(logical(GateArray(use_ind,:))) = 1;

for i = 1:size(GateArray,1)
    weights(logical(GateArray(i,:))) = exp(-abs(i-use_ind));
end

