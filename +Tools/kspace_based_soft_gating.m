function weights = kspace_based_soft_gating(GateArray,use_ind)

weights(GateArray==use_ind) = 1;

for i = 1:length(GateArray)
    weights((GateArray==i)) = exp(-2*abs(i-use_ind));
end
