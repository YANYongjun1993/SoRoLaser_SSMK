load figure8_lowDim_trajectory.mat
ref = struct;
ref.y = lowDim';
ref.t = 0.01:0.01:size(lowDim,2)*0.01;
save("figure8_lowDim_traj.mat","ref")