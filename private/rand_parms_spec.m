function [parms] = rand_parms_spec()
parms.D       =  3;
parms.K       =  3;
parms.alpha   =  .6;            % balancing parameter: 1 perfectly balanced
parms.r       =  0.1;            % outlier ratio
parms.sigma   =  0;           % stdv of Gaussian noise
return