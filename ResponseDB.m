%> @file  ResponseDB.m
%> @brief Response function database.
% ==============================================================================
%> @brief Read and write response data.
%
%> This class reads and writes response data to/from an HDF5 file.  It 
%> populates response operators as functions of keff for each node type and
%> provides vectorized interpolation.
% ==============================================================================
classdef ResponseDB < handle
    
    properties (Constant)

    end

    properties (Access = private)
        d_input
        d_k_vector
        d_file
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
            this.d_k_vector = get(input, 'rf_k_vector');
            
            create(this, '/keffs', [length(this.d_k_vector) 1]);
            h5writeatt(file, '/', 'dimension',          get(input, 'dimension'));
            h5writeatt(file, '/', 'solver',             'matlab-transport-pack');
            h5writeatt(file, '/', 'number_groups',      get(input, 'number_groups'));
            h5writeatt(file, '/', 'max_order_space',    get(input, 'rf_max_order_space'));
            h5writeatt(file, '/', 'max_order_azimuth',  get(input, 'rf_max_order_azimuth'));
            h5writeatt(file, '/', 'max_order_polar',    get(input, 'rf_max_order_polar'));
            h5writeatt(file, '/', 'number_keff',        length(this.d_k_vector));
            h5write(file, '/keffs', this.d_k_vector);
            this.d_file = file;
        end
        
        % ======================================================================
        %> @brief Write responses
        %>
        %> 
        %>
        %> @param R     Cell array of response blocks for all keffs.
        % ======================================================================        
        function this = write_response(this, node_index, node_description, ...
                R, F, A, L)
            
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
        