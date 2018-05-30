function Fit_Example_Data()

%% ======================== DATA INPUTS ========================

% First, read in the parentage matrix. Each row corresponds to a source reef containing
% sampled parents. Each column corresponds to a destination reef containing sampled
% juveniles. The elements in each cell indicate the number of assigned juveniles on the
% COLUMN reef whose parents were on the ROW reef. The final row contains the unassigned
% juveniles sampled on each COLUMN reef.
Assignments = csvread('Example_parentage_matrix.csv');

% Second, read in the proportion of adults on each reef that were sampled. This vector
% should have the number of rows equal to the number of columns in the "Obs" matrix.
% That is - the number of sampled reefs.
Adult_sample_proportions = csvread('Example_proportion_sampled.csv');

% Third, read in a distance matrix showing the distances between all
Distances = csvread('Example_distance_matrix.csv');

% Fourth, read in a distance matrix showing the distances between all
Reef_sizes = csvread('Example_patch_size_matrix.csv');

% Fifth, read in a list of the reefs that were sampled
Sampled_reefs = csvread('Example_sampled_patch_list.csv');

% Sixth, read in the locations of each of the reefs
Centroids = csvread('Example_patch_location.csv');

%% =============================================================

% How many reefs were sampled
Num_sampled_reefs = size(Assignments,2);
Num_reefs = size(Centroids,1);

% Define the generalised Gaussian functions that we're fitting here
F  = @(x,k,theta)    exp(k).*theta.*exp(-(exp(k)*x).^theta)/gamma(1/theta);
FM = @(x,k,theta) x.*exp(k).*theta.*exp(-(exp(k)*x).^theta)/gamma(1/theta);

% Run through the list of potential kernels and fit each one to the data
Theta_list = [1 2 3 0.5];
LowerBound = -10;
UpperBound = 10;
for th = 1:length(Theta_list);
   [Best_k(th),LL_k(th)] = fminbnd(@Kernel_Fitting_Function,LowerBound,UpperBound,[],... % These are the search input parameters
      Assignments,Distances,Reef_sizes,Sampled_reefs,Adult_sample_proportions,F,Theta_list(th)); % These are the extra parameters needed by the function
end
if min(Best_k) < 1.01.*LowerBound
   disp('Optimisation is choosing lower bound - reduce value to avoid error')
elseif max(Best_k) > 0.99.*UpperBound
   disp('Optimisation is choosing upper bound - increase value to avoid error')
end

% Identify the best fit from the candidate kernels
[~,Best_kernel] = min(LL_k);
Best_fit_k = Best_k(Best_kernel);

% Create bootstrap confidence bounds around the best estimate by re-sampling at the patch scale
for b = 1:100
   
   % First resample the assignment data, as well as the proportion of each reef sampled
   Bootstrap_resample = randsample(Num_sampled_reefs,Num_sampled_reefs,1);
   Sampled_reefs_b = Sampled_reefs(Bootstrap_resample);
   Adult_sample_proportions_b = Adult_sample_proportions(Bootstrap_resample);
   Assignments_b = [Assignments(Bootstrap_resample,Bootstrap_resample); Assignments(end,:)];
   
   % Re-fit the kernel of the correct shape to this data
   [Bootstrap_k(b),LL_b] = fminbnd(@Kernel_Fitting_Function,LowerBound,UpperBound,[],... % These are the search input parameters
      Assignments_b,Distances,Reef_sizes,Sampled_reefs_b,Adult_sample_proportions_b,F,Theta_list(Best_kernel)); % These are the extra parameters needed by the function
end
Confidence_bounds = quantile(Bootstrap_k,[0.975 0.025]);

% Calculate the mean dispersal distance by integrating the best fit kernel, 
% and the confidence intervals from the quantiles of the bootstrap re-sampled fits.
MDD   =  integral(@(x)FM(x,Best_fit_k,Theta_list(Best_kernel)),0,500);
CI_DD = [integral(@(x)FM(x,Confidence_bounds(1),Theta_list(Best_kernel)),0,500) ...
         integral(@(x)FM(x,Confidence_bounds(2),Theta_list(Best_kernel)),0,500)];

disp(['Mean dispersal distance is ' num2str(MDD,3) ' km, with 95% confidence intervals of [' num2str(CI_DD(1),3) ', ' num2str(CI_DD(2),3) '] km'])

if Best_fit_k == LowerBound
   disp('ERROR: Reduce the lower bound passed to function FMINBND')
elseif Best_fit_k == UpperBound
   disp('ERROR: Increase the upper bound passed to function FMINBND')
end

save OUTPUTS








