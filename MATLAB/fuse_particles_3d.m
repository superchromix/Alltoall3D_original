% fuse_particles fuses a set of particles by translation and rotation
% transformations of their 3D coordinates.
%
% Parameters:
%   Input:
%       n_particles: number of particles
%       n_localizations_per_particle: number of localizations for each
%                                     particle
%       coordinates_x: x coordinates of each localization
%       coordinates_y: y coordinates of each localization
%       coordinates_z: z coordinates of each localization
%       precision_xy: localization uncertainties in xy for each localization
%       precision_z: localization uncertainties in z for each localization
%       gauss_transform_scale: ??
%       channel_ids: the channel ids of each localization (optional,
%                    default values are 0)
%       averaging_channel_id: channel_id of the channel the fusion is based
%                             on (optional, default value is 0)
%       n_iterations_all2all: number of iterations of the alltoall routine
%                             (optional, dafault value is 1)
%       n_iterations_one2all: number of iterations of the onetoall routine
%                             (optional, dafault value is 10)
%       symmetry_order:       the symmetry order in the xy plane (optional)
%       outlier_threshold: threshold used for the outlier removal, low
%                          value -> strong removal (optional, default value
%                          is 1)
%   Output:
%       transformed_coordinates_x: transformed x coordinates of each
%                                  localization
%       transformed_coordinates_y: transformed y coordinates of each
%                                  localization
%       transformed_coordinates_z: transformed z coordinates of each
%                                  localization
%       transformation_parameters: final transformation parameters for each
%                                  particle

function [transformed_coordinates_x, transformed_coordinates_y, transformed_coordinates_z, transformation_parameters]...
    = fuse_particles_3d(...
        n_particles,...
        n_localizations_per_particle,...
        coordinates_x,...
        coordinates_y,...
        coordinates_z,...
        precision_xy,...
        precision_z,...
        gauss_transform_scale,...
        channel_ids,...
        averaging_channel_id,...
        n_iterations_all2all,...
        n_iterations_one2all,...
        symmetry_order,...
        outlier_threshold)

%% GPU flags
USE_GPU_EXPDIST = true;
USE_GPU_GAUSSTRANSFORM = true;

if gpuDeviceCount <= 0
    USE_GPU_GAUSSTRANSFORM = false;
    USE_GPU_EXPDIST = false;
end

if USE_GPU_GAUSSTRANSFORM 
    if ~exist('mex_gausstransform','file')
        USE_GPU_GAUSSTRANSFORM = false;
    end
end
    
if USE_GPU_EXPDIST 
    if ~exist('mex_expdist','file')
        USE_GPU_EXPDIST = false;
    end
end

if ~USE_GPU_GAUSSTRANSFORM
    if ~exist('mex_gausstransform_cpu','file')
        message = 'No compiled modules found for GaussTransform.\n';
        message_id = 'MATLAB:MEXNotFound';
        error (message_id, message);
    end
    fprintf('running CPU version of gausstransform\n')
else
    fprintf('running GPU version of gausstransform\n')
end

if ~USE_GPU_EXPDIST 
    if ~exist('mex_expdist_cpu','file')
        message = 'No compiled modules found for ExpDist.\n';
        message_id = 'MATLAB:MEXNotFound';
        error (message_id, message);
    end
    fprintf('running CPU version of expdist\n')
else
    fprintf('running GPU version of expdist\n')
end

%% setting indicies of the first localization of each particle
particle_beginnings = ones(n_particles,1);
particle_endings(n_particles) = numel(coordinates_x);
for i = 2:n_particles
    particle_beginnings(i) = particle_beginnings(i-1) + n_localizations_per_particle(i-1);
    particle_endings(i-1) = particle_beginnings(i)-1;
end

%% starting parallel pool
pp = gcp;
if ~(pp.Connected)
    parpool();
end

%% setting channel filter
channel_filter = channel_ids == averaging_channel_id;

%% performing the all2all registration
pprint('all2all registration ',45)
t = tic;
coordinates = [coordinates_x, coordinates_y, coordinates_z];
precision = [precision_xy, precision_z];
all2all_matrix = cell(n_particles-1,n_particles);
particle_begin_i = 1;
coordinates_j = cell(n_particles,1);
precision_j = cell(n_particles,1);
for i=1:n_particles-1
    
    indices_i = particle_beginnings(i):particle_endings(i);
    indices_i = indices_i(channel_filter(indices_i));
    coordinates_i = coordinates(indices_i,:);
    precision_i = precision(indices_i,:);
    
    for j=i+1:n_particles
        indices_j = particle_beginnings(j):particle_endings(j);
        indices_j = indices_j(channel_filter(indices_j));
        coordinates_j{j} = coordinates(indices_j,:);
        precision_j{j} = precision(indices_j,:);
    end
    
    parfor j=i+1:n_particles
        all2all_matrix{i,j}.parameters = all2all3Dn(coordinates_i, coordinates_j{j}, precision_i, precision_j{j}, n_iterations_all2all, gauss_transform_scale, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
        all2all_matrix{i,j}.ids = [i; j];
    end
    
    progress_bar(n_particles-1,i);
    
    particle_begin_i = particle_begin_i + n_localizations_per_particle(i);
end
fprintf([' ' num2str(toc(t)) ' s\n']);

%% translating quaternion to matrix representation (SE3)
k = 1;
RR = zeros(4,4,n_particles*(n_particles-1)/2);
I = zeros(2,n_particles*(n_particles-1)/2);
for i=1:n_particles-1
    for j=i+1:n_particles
        q = [all2all_matrix{i,j}.parameters(4) all2all_matrix{i,j}.parameters(1) ...
             all2all_matrix{i,j}.parameters(2) all2all_matrix{i,j}.parameters(3)];
         
        % RR holds the registration parameters of size 4x4xn_particles(n_particles-1)/2 
        RR(1:3,1:3,k) = q2R(q);
        RR(1:3,4,k) = [all2all_matrix{i,j}.parameters(5); all2all_matrix{i,j}.parameters(6); all2all_matrix{i,j}.parameters(7)];
        RR(4,4,k) = 1;
        RR(4,1:3) = 0;
        
        % I holds the connectivity of pairs of size 2xn_particles(n_particles-1)/2
        I(:,k)=all2all_matrix{i,j}.ids;
        k=k+1;        
    end
end

%% averaging transformation parameters
pprint('averaging transformation parameters ',45);
t = tic;
Mest = MeanSE3Graph(RR, I);
progress_bar(1,1);
fprintf([' ' num2str(toc(t)) ' s\n']);

if outlier_threshold < pi
    %% remove outliers
    pprint('removing outliers ',45);
    t = tic;
    [RR, I] = remove_outliers(RR, I, Mest, outlier_threshold);
    progress_bar(1,1);
    fprintf([' ' num2str(toc(t)) ' s\n']);

    %% repeating averaging transformation parameters
    pprint('averaging transformation parameters ',45);
    t = tic;
    Mest = MeanSE3Graph(RR, I);
    progress_bar(1,1);
    fprintf([' ' num2str(toc(t)) ' s\n']);
end

%% applying the absolute registration parameters to particles
pprint('coordinate transformation ',45);
t = tic;
transformed_coordinates = zeros(sum(n_localizations_per_particle),3,n_iterations_one2all+1);
transformed_particles = cell(1,n_particles);
for i=1:n_particles  
    
    indices = particle_beginnings(i):particle_endings(i);
    indices = indices(channel_filter(indices));
    
    estA = eye(4);
    estA(1:3,1:3) = Mest(1:3,1:3,i); 
    estA(4,:) = Mest(:,4,i)';
    estTform = affine3d(estA);
    
    % transform coordinates
    coordinates_ptc = pointCloud(coordinates(indices,:));
    transformed_coordinates_ptc = pctransform2(coordinates_ptc, invert(estTform));
    
    % copy transformed coordinates
    transformed_particles{i}.points = transformed_coordinates_ptc.Location;
    transformed_particles{i}.sigma = [precision_xy(indices), precision_z(indices)];
    transformed_coordinates(indices,:,1) = transformed_coordinates_ptc.Location;
end
progress_bar(1,1);
fprintf([' ' num2str(toc(t)) ' s\n']);

%% performing the one2all registration
pprint('one2all registration ',45);
t = tic;
tc = one2all3D(transformed_particles, n_iterations_one2all, [], '.', transformed_coordinates(channel_filter,:,1), gauss_transform_scale, symmetry_order, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
transformed_coordinates(channel_filter,:,:) = reshape(cell2mat(tc),[],3,n_iterations_one2all+1);
fprintf([' ' num2str(toc(t)) ' s\n']);

%% calculationg final transformation parameters
transformation_parameters = zeros(4,4,n_particles,n_iterations_one2all+1);
for j = 1:n_iterations_one2all+1
    for i = 1:n_particles

        indices = particle_beginnings(i):particle_endings(i);
        indices = indices(channel_filter(indices));

        transformation_parameters(:,:,i,j) = get_final_transform_params(...
            transformed_coordinates(indices,:,j), coordinates(indices,:));
    end
end

%% transforming remaining channels
for iter = 1:n_iterations_one2all+1
    for i = 1:n_particles
        tp.rot = transformation_parameters(1:3,1:3,i,iter);
        tp.shift = transformation_parameters(4,1:3,i,iter);

        indices = particle_beginnings(i):particle_endings(i);
        indices = indices(~channel_filter(indices));

        transformed_coordinates(indices,:,iter) =  coordinates(indices,:) * tp.rot + tp.shift;
    end
end

%% setting output arguments holding transformed coordinates
transformed_coordinates_x = squeeze(transformed_coordinates(:,1,:));
transformed_coordinates_y = squeeze(transformed_coordinates(:,2,:));
transformed_coordinates_z = squeeze(transformed_coordinates(:,3,:));

end
