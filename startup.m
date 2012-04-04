%> @file  startup.m
%> @brief Set paths and other things.

% User configuration options.
use_petsc = 1;
use_slepc = 1;
mtp_dir   = '/home/robertsj/Research/matlab_transport_pack/source';
% turn off threading for testing
maxNumCompThreads(1);


% Include source directories.
path(path, mtp_dir)
path(path, [mtp_dir,'/utilities'])
path(path, [mtp_dir,'/angle'])
path(path, [mtp_dir,'/boundary'])
path(path, [mtp_dir,'/diffusion'])
path(path, [mtp_dir,'/erme'])
path(path, [mtp_dir,'/geometry'])
path(path, [mtp_dir,'/material'])
path(path, [mtp_dir,'/solvers'])
path(path, [mtp_dir,'/transport'])
path(path, [mtp_dir,'/tests'])

% PETSc and SLEPc
if use_petsc && ~exist('PetscInitialize','file')
  PETSC_DIR = getenv('PETSC_DIR');
  if isempty(PETSC_DIR) 
    error('Must set environment variable PETSC_DIR or add the appropriate dir to Matlab path')
  end
  disp(['Using PETSC_DIR=', PETSC_DIR])
  path(path,[PETSC_DIR '/bin/matlab/classes'])
end
if use_slepc && ~exist('SlepcInitialize','file')
  SLEPC_DIR = getenv('SLEPC_DIR');
  if isempty(SLEPC_DIR) 
    error('Must set environment variable SLEPC_DIR or add the appropriate dir to Matlab path')
  end
  disp(['Using SLEPC_DIR=', SLEPC_DIR])
  path(path,[SLEPC_DIR '/bin/matlab/classes'])
end
