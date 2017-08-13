function LL_k = Kernel_Fitting_Function(k,Obs,Dist,Area,ID,SS_A,Kernel,POWER)
 
% ----- NAME -----
% 
% Kernel_Fitting_Function 
% 
% Updated 2 June 2017 by Michael Bode. Email: michael.bode@jcu.edu.au
% 
% ----- DESCRIPTION -----
% 
% Calculates the likelihood of a particular observed parentage dataset, given 
% a kernel functional form, and a particular parameterisation. The function 
% implements the likelihood formula given in Eq. 5 in Bode et al. (Submitted), 
% "Estimating dispersal kernels using genetic parentage data", submitted to 
% Methods in Ecology and Evolution
% 
% ----- INPUTS -----
% 
% k
%    A vector of parameters for the inline function (the dispersal kernel)
%    NB: this input is only the parameter being optimised. 
% 
% Obs
%    The observed parentage data. This is the matrix M described in the manuscript.
%    
% Dist
%    This is a PxP square matrix containing the distance between all patches 
%    in the analysis region. (P is the number of patches.)
%    
% Area
%    This is a Px1 vector containing the area of all patches in the analysis region
% 
% ID
%    This is a vector containing the identity (i.e., the index number) of the sampled 
%    patches. It is the set S_J described in the manuscript.
% 
% SS_A
%    This is a Px1 vector containing the proportion of adults sampled on each patch
%    in the analysis region. It contains the variables \pi_i in the manuscript.
% 
% Kernel
%    This is a matlab inline function that describes the dispersal kernel, e.g.:
%       F  = @(x,k,theta)    exp(k).*theta.*exp(-(exp(k)*x).^theta)/gamma(1/theta);
%    This is equivalent to Eq. 1 in the manuscript.
% 
% POWER
%    This is the shape function in a generalised Gaussian dispersal kernel, if needed.
%    It is \theta in Eq. 1 in the manuscript. 
% 
% ----- OUTPUTS -----
% 
% LL_k
%    This is the likelihood of the dataset OBS, given parameter P.
%    
 
% Grab some parameters from the inputs
NumReefs = length(Dist);
 
% For each reef, create the expected proportions between every reef pair
Proportions = Kernel(Dist, k, POWER);
 
% Inflate the values by the total area of the source reef
% (larger reefs send out more larvae, all else being equal)
Settlers = Proportions.*repmat(Area,1,NumReefs);
 
% Sum the total settlers on each reef
AllSettlers = sum(Settlers);
 
% How many reefs were sampled?
NumSampledReefs = length(ID);
 
% Parentage is assigned when both the source and sink reef have been sampled.
AssignedSettlers = zeros(NumSampledReefs+1,NumSampledReefs);
for i = 1:NumSampledReefs
   This_SS_A = SS_A(i);
   for j = 1:NumSampledReefs
      SettlersFromAssignedReefs = Settlers(ID(i),ID(j));
      
      % Not all settlers from assigned reefs will be assigned, because not all adults were sampled
      AssignedSettlers(i,j) = SettlersFromAssignedReefs.*(This_SS_A.^2 + 2.*This_SS_A.*(1 - This_SS_A));
      AssignedSettlers(NumSampledReefs+1,j) = AssignedSettlers(NumSampledReefs+1,j) ...
         + SettlersFromAssignedReefs.*(1-This_SS_A).^2;
   end
end
 
% Parentage is unknown (the last row of the matrix) when the settlers
% come from an unsampled reef but disperse to a sampled reef.
Unsampled = setdiff(1:NumReefs,ID);
for j = 1:NumSampledReefs
   AssignedSettlers(NumSampledReefs+1,j) = AssignedSettlers(NumSampledReefs+1,j) ...
      + sum(Settlers(Unsampled,ID(j)));
end
 
% Normalise values into multinomial probabilities
PredictedProportions = AssignedSettlers./repmat(sum(AssignedSettlers),NumSampledReefs+1,1);
 
% Loglikelihoods can't handle zeros. Make them very small instead.
PredictedProportions(PredictedProportions==0) = 1e-12;
 
% Re-normalise into probabilities
PredictedProportions = PredictedProportions./repmat(sum(PredictedProportions),size(PredictedProportions,1),1);
 
LL_k = 0; % Initialise
for i = 1:NumSampledReefs 
   % Go through reefs one-by-one: What is the LL of each reef's observations?
   
   ObsVector = Obs(:,i);
   ProbVector = PredictedProportions(:,i);
   
   % Calculate the log likelihood for this particular reef, given the sample that's been observed
   LL_k = LL_k + sum(ObsVector.*log(ProbVector));
   
end
 
% Invert the likelihood if using a matlab routine that minimises (e.g., FMINCON, FMINSEARCH, FMINBND)
LL_k = -LL_k;
 
