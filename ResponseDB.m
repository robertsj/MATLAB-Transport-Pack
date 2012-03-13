%> @file  ResponseDB.m
%> @brief Response function database.
% ==============================================================================
%> @brief Read and write response data.
%
%> This class reads and writes response data to/from an HDF5 file.  It 
%> populates response operators as functions of keff for each node type and
%> provides vectorized interpolation.
%>
%> @note requires "hdf5tools" from 
%>       http://www.mathworks.com/matlabcentral/fileexchange/17172
% ==============================================================================
classdef ResponseDB < handle
    
    properties (Constant)

    end

    properties (Access = private)
        d_input
        d_k_vector
        d_file
        d_file_id
        %
        d_number_nodes
        d_number_keffs
        d_number_groups
        d_max_order_space
        d_max_order_azimuth
        d_max_order_polar
        %> Responses for all elements and all keffs
        d_R
        d_F
        d_A
        d_L
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief  Class constructor
        %> @param  Input Input database.
        %> @return Instance of the ResponseDB class.
        % ======================================================================
        function this = ResponseDB(input) 
            this.d_input = input;
            file = get(input, 'rf_db_name');
            if (~file)
                file = 'default.h5';
            end
            this.d_file = file;
        end
        
        % ======================================================================
        %> @brief  Initialize HDF5 file for writing.
        % ======================================================================
        function this = initialize_write(this)
            
            % Open file id.  Overwrite the file if it exists.
            this.d_file_id = h5filecreate(this.d_file, 'truncate', true);
            H5F.close(this.d_file_id);
            
            this.d_k_vector = get(this.d_input, 'rf_k_vector');
            
            create(this, '/keffs', [length(this.d_k_vector) 1]);
            h5writeatt(this.d_file, '/', 'dimension',          get(this.d_input, 'dimension'));
            h5writeatt(this.d_file, '/', 'solver',             'matlab-transport-pack');
            h5writeatt(this.d_file, '/', 'number_groups',      get(this.d_input, 'number_groups'));
            h5writeatt(this.d_file, '/', 'max_order_space',    get(this.d_input, 'rf_max_order_space'));
            h5writeatt(this.d_file, '/', 'max_order_azimuth',  get(this.d_input, 'rf_max_order_azimuth'));
            h5writeatt(this.d_file, '/', 'max_order_polar',    get(this.d_input, 'rf_max_order_polar'));
            h5writeatt(this.d_file, '/', 'number_keffs',       length(this.d_k_vector));
            h5writeatt(this.d_file, '/', 'number_nodes',       get(this.d_input, 'rf_number_nodes'));
            h5write(this.d_file, '/keffs', this.d_k_vector);

        end
        
        % ======================================================================
        %> @brief Write responses
        %
        %> @param R     Cell array of response blocks for all keffs.
        % ======================================================================        
        function this = write_response(this, node_index, node_description, ...
                R, F, A, L)
          	if length(this.d_k_vector) == 1
                warning('user:input','Did you mean to produce a db with one k?')
            end
            % convert node index to string
            nidx = ['/n',num2str(node_index)];
            size_R = size(R(:, :, 1));
            size_F = size(F(:, 1));
            size_A = size(A(:, 1));
            size_L = size(L(:, :, 1));
            for k_index = 1:length(this.d_k_vector)
                location = [nidx,'/k',num2str(k_index)];
                create(this, [location,'/R'], size_R);
                create(this, [location,'/F'], size_F);
                create(this, [location,'/A'], size_A);
                create(this, [location,'/L'], size_L);
                h5write(this.d_file, [location,'/R'], R(:, :, k_index));
                h5write(this.d_file, [location,'/F'], F(:, k_index));
                h5write(this.d_file, [location,'/A'], A(:, k_index));
                h5write(this.d_file, [location,'/L'], L(:, :, k_index));
            end
            h5writeatt(this.d_file, nidx, 'description', node_description);
            h5writeatt(this.d_file, nidx, 'width_x', 1.0);
            h5writeatt(this.d_file, nidx, 'width_y', 1.0);
        end
        
        % ======================================================================
        %> @brief Read responses.
        % ======================================================================        
        function this = read_response(this)
            
            % Read attributes
            this.d_number_nodes      = h5readatt(this.d_file, '/', 'number_nodes');
            if this.d_number_nodes == 0
                error('Read zero nodes!')
            end
            this.d_number_keffs      = h5readatt(this.d_file, '/', 'number_keffs');
            if this.d_number_keffs == 0
                error('Read zeros keffs!')
            end
            this.d_number_groups     = h5readatt(this.d_file, '/', 'number_groups');
            this.d_max_order_space   = h5readatt(this.d_file, '/', 'max_order_space');
            this.d_max_order_azimuth = h5readatt(this.d_file, '/', 'max_order_azimuth');
            this.d_max_order_polar   = h5readatt(this.d_file, '/', 'max_order_polar');
            
            % Read vector of keffs
            this.d_k_vector = h5read(this.d_file, '/keffs');            
            
            len_R = 4 * this.d_number_groups * (1+this.d_max_order_space) * ...
                (1+this.d_max_order_azimuth) * (1+this.d_max_order_polar);
            this.d_R = zeros(len_R, len_R, this.d_number_keffs, ...
                this.d_number_nodes);
            this.d_F = zeros(len_R, this.d_number_keffs, this.d_number_nodes);
            this.d_A = zeros(len_R, this.d_number_keffs, this.d_number_nodes);
            this.d_L = zeros(len_R, 4, this.d_number_keffs, ...
                this.d_number_nodes);            
            
            % convert node index to string
            for ni = 1:this.d_number_nodes
                nidx = ['/n',num2str(ni)];
                for ki = 1:this.d_number_keffs
                    kidx = ['/k',num2str(ki)];
                    loc = [nidx, kidx];
                    this.d_R(:, :, ki, ni) = h5read(this.d_file, [loc, '/R']);
                    this.d_F(:,    ki, ni) = h5read(this.d_file, [loc, '/F']);
                    this.d_A(:,    ki, ni) = h5read(this.d_file, [loc, '/A']);
                    this.d_L(:, :, ki, ni) = h5read(this.d_file, [loc, '/L']);
                end
            end

        end        
        
        % ======================================================================
        %> @brief Return all responses (for all keffs and nodes).
        % ======================================================================            
        function [R, F, A, L] = get_all_responses(this)
           R = this.d_R;
           F = this.d_F;
           A = this.d_A;
           L = this.d_L;
        end
        
        % ======================================================================
        %> @brief Return responses for all nodes and an interpolated keff.
        %> @param   keff            Value of keff to interpolate responses.
        %> @param   interp_method   Method of interpolation (default cubic).
        %> @return                  Responses    
        % ======================================================================    
        function [R, F, A, L] = get_responses(this, keff, interp_method)
            if nargin == 2
                interp_method = 'spline'; % linear is fast, but no extrapolation
            end
            len_R = 4 * this.d_number_groups * (1+this.d_max_order_space) * ...
                (1+this.d_max_order_azimuth) * (1+this.d_max_order_polar);
            nk = this.d_number_keffs;    
            % Interpolated values
            R = zeros(len_R, len_R, this.d_number_nodes);
            F = zeros(len_R,        this.d_number_nodes);
            A = zeros(len_R,        this.d_number_nodes);
            L = zeros(len_R,     4, this.d_number_nodes);
            for ni = 1:this.d_number_nodes;
                tmpR = reshape(this.d_R(:, :, :, ni), len_R^2, nk);
                R(:, :, ni) = reshape( ...
                    interp1(this.d_k_vector, tmpR', keff, interp_method), ...
                    len_R, len_R);
                tmpF(:, :) = this.d_F(:, :, ni)';
                F(:, ni) = interp1(this.d_k_vector, tmpF, keff, interp_method)';
                tmpA(:, :) = this.d_A(:, :, ni)';
                A(:, ni) = interp1(this.d_k_vector, tmpA, keff, interp_method)';      
                tmpL = reshape(this.d_L(:, :, :, ni), len_R*4, nk);
                L(:, :, ni) = reshape( ...
                    interp1(this.d_k_vector, tmpL', keff, interp_method), ...
                    len_R, 4);
            end
            
        end        
        
    end
    
    methods (Access = private)
       
        % ======================================================================
        %> @brief Create a dataset.
        %>
        %> This does nothing if the data is already there.
        %>
        %> @param   location    Location with respect to root.
        %> @param   datasize    Size of data array
        % ======================================================================        
        function this = create(this, location, datasize, overwrite)
            if nargin == 3
                overwrite = 1;
            end
            try
                h5create(this.d_file, location, datasize);
            catch me
                disp(me.message);
                idexists = 'MATLAB:imagesci:h5create:datasetAlreadyExists'; 
                if strcmp(me.identifier, idexists) && ~overwrite
                    error(['Data set ', location, ' in file ', this.d_file, ...
                          ' already exists and you said not to overwrite!'])
                elseif ~strcmp(me.identifier, idexists)
                    error(['Unanticipated error creating dataset: ', ...
                           location, ' in file ', this.d_file]);
                end
            end
        end
        
    end
end
        