%> @file  ResponseDB.m
%> @brief Response function database.
% ==============================================================================
%> @brief Read and write response data.
%
%> This class reads and writes response data to/from an HDF5 file.  It 
%> populates response operators as functions of keff for each node type and
%> provides vectorized interpolation.
%>
%> Relevant input database entries include thos of @ref ResponseServerBase
%> and the following
%>   - rf_db_name
%>   - rf_keff_vector
%>   - rf_max_order_space
%>   - rf_max_order_azimuth
%>   - rf_max_order_polar
%>   - rf_interp_method
%>
%> @note Requires "hdf5tools" from 
%>       http://www.mathworks.com/matlabcentral/fileexchange/17172
%>       Only one function is used now, but more may be used later.
%>
%> @note Most stuff is specific to 2-d, but I'm trying to build in a more 
%>       general approach for future 3-d testing.  
%>
%> @todo I only have basic R, F, A, and L operators.  Eventually, power edit
%>       regions will be necessary.  However, those are almost a post-
%>       processing thing, as they won't affect the solution and need not be
%>       in memory before hand.
%>
% ==============================================================================
classdef ResponseDB < ResponseServerBase

    properties (Access = protected)
        
        %> DB file name and id
        d_file
        d_file_id
        
        %> Cell array of vectors of node keffs.
        d_keff_vector
        
        %> Array of the number of keffs for each node.
        d_number_keffs
        
        %> Cell array of string descriptions fo reach node.
        d_node_descriptions
        
        %> Interpolation method
        d_interp_method
        
        %> Responses for all elements and all keffs.  These are cell arrays
        %> of each set of responses per node.
        d_R_all
        d_F_all
        d_A_all
        d_L_all
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief  Class constructor
        %> @param  Input Input database.
        %> @return Instance of the ResponseDB class.
        % ======================================================================
        function this = ResponseDB(input) 
            
            % Call base class
            this = this@ResponseServerBase(input);
            
            % HDF5 file
            file = get(input, 'rf_db_name');
            if (~file)
                file = 'default.h5';
            end
            this.d_file = file;
            % Default the interpolation to cubic, which seems to be a lot better
            % than linear and much, much faster than cubic (even though the
            % documentation says they should cost the same).
            this.d_interp_method = get(input, 'rf_interp_method');
            if ~this.d_interp_method 
                this.d_interp_method = 'spline';
            end
        end
        
        % ======================================================================
        %> @brief  Initialize HDF5 file for writing.
        % ======================================================================
        function this = initialize_write(this)
            
            % Open file id.  Overwrite the file if it exists.
            this.d_file_id = h5filecreate(this.d_file, 'truncate', true);
            H5F.close(this.d_file_id);
            
            keffs  = get(this.d_input, 'rf_keff_vector');
            this.d_number_nodes = get(this.d_input, 'rf_number_nodes');
            DBC.Require('length(keffs)==this.d_number_nodes');
            this.d_keff_vector = keffs;
            
            db_description = get(this.d_input, 'rf_db_description');
            if ~db_description
                db_description = 'just another database.';
            end
            
            h5writeatt(this.d_file, '/', 'db_description',          ...
                get(this.d_input, 'rf_db_description'));
            h5writeatt(this.d_file, '/', 'dimension',          ...
                get(this.d_input, 'dimension'));
            h5writeatt(this.d_file, '/', 'solver',             ...
                'matlab-transport-pack');
            h5writeatt(this.d_file, '/', 'number_groups',      ...
                get(this.d_input, 'number_groups'));
            h5writeatt(this.d_file, '/', 'max_order_space',    ...
                get(this.d_input, 'rf_max_order_space'));
            h5writeatt(this.d_file, '/', 'max_order_azimuth',  ...
                get(this.d_input, 'rf_max_order_azimuth'));
            h5writeatt(this.d_file, '/', 'max_order_polar',    ...
                get(this.d_input, 'rf_max_order_polar'));
            h5writeatt(this.d_file, '/', 'number_keffs',       ...
                length(this.d_keff_vector));
            h5writeatt(this.d_file, '/', 'number_nodes',       ...
                this.d_number_nodes);
           

        end
        
        % ======================================================================
        %> @brief Write responses
        %
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
            
            % Get this node's number of keffs and write them.
            keffs  = this.d_keff_vector{node_index};
            size_k = length(this.d_keff_vector{node_index});
            create(this, [nidx,'/keffs'], size_k);
            h5write(this.d_file, [nidx,'/keffs'], keffs);
            h5writeatt(this.d_file, nidx, 'number_keffs', size_k);
            
            % Write the responses for each keff.
            for k_index = 1:size_k
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
            
            % Other attributes.
            h5writeatt(this.d_file, nidx, 'description', node_description);
            h5writeatt(this.d_file, nidx, 'width_x', 1.0);
            h5writeatt(this.d_file, nidx, 'width_y', 1.0);

        end
        
        % ======================================================================
        %> @brief Read responses.
        % ======================================================================        
        function this = read_response(this)
            
            f = this.d_file;
            
            % Read attributes
            this.d_number_nodes = h5readatt(f, '/', 'number_nodes');
            if this.d_number_nodes == 0
                error('Read zero nodes!')
            end

            this.d_dimension         = h5readatt(f, '/', 'dimension');
            this.d_number_groups     = h5readatt(f, '/', 'number_groups');
            this.d_max_order_space   = h5readatt(f, '/', 'max_order_space');
            this.d_max_order_azimuth = h5readatt(f, '/', 'max_order_azimuth');
            this.d_max_order_polar   = h5readatt(f, '/', 'max_order_polar');
            
            % Check whether these maximum orders are less than the requested
            % orders.  If they are, reset the requested orders to match.
            if this.d_max_order_space < this.d_order_space
                error('user:input', 'Requested spatial order > max');
                %this.d_order_space = this.d_max_order_space;
            end
            if this.d_max_order_azimuth < this.d_order_azimuth
                error('user:input', 'Requested azimuth order > max');
                %this.d_order_azimuth = this.d_max_order_azimuth;
            end
            if this.d_max_order_polar < this.d_order_polar
                error('user:input', 'Requested polar order > max');
                %this.d_order_polar = this.d_max_order_polar;
            end            

            % Block size
            len_R = 4 * this.d_number_groups * (1+this.d_max_order_space) * ...
                (1+this.d_max_order_azimuth) * (1+this.d_max_order_polar);
            
            this.d_R_all = cell(this.d_number_nodes, 1);
            this.d_F_all = cell(this.d_number_nodes, 1);
            this.d_A_all = cell(this.d_number_nodes, 1);
            this.d_L_all = cell(this.d_number_nodes, 1);
            this.d_keff_vector = cell(this.d_number_nodes, 1);
            this.d_number_keffs = zeros(this.d_number_nodes, 1);
            for ni = 1:this.d_number_nodes
                
                nidx = ['/n',num2str(ni)];
                
                % Read node description.
                this.d_node_descriptions{ni} = h5readatt(f, nidx, 'description');
                
                % Read keffs.
                this.d_number_keffs(ni) = h5readatt(f, nidx, 'number_keffs');
                this.d_keff_vector{ni} = h5read(f, [nidx, '/keffs']);
                
                this.d_R_all{ni} = zeros(len_R, len_R, this.d_number_keffs(ni));
                this.d_F_all{ni} = zeros(len_R,        this.d_number_keffs(ni));
                this.d_A_all{ni} = zeros(len_R,        this.d_number_keffs(ni));
                this.d_L_all{ni} = zeros(len_R,     4, this.d_number_keffs(ni));
                
                for ki = 1:this.d_number_keffs(ni)
                    kidx = ['/k',num2str(ki)];
                    loc = [nidx, kidx];
                    this.d_R_all{ni}(:, :, ki) = h5read(f, [loc, '/R']);
                    this.d_F_all{ni}(:,    ki) = h5read(f, [loc, '/F']);
                    this.d_A_all{ni}(:,    ki) = h5read(f, [loc, '/A']);
                    this.d_L_all{ni}(:, :, ki) = h5read(f, [loc, '/L']);
                end
            end

        end        
        
        % ======================================================================
        %> @brief Return all responses (for all keffs and nodes).
        % ======================================================================            
        function [R, F, A, L] = get_all_responses(this)
           R = this.d_R_all;
           F = this.d_F_all;
           A = this.d_A_all;
           L = this.d_L_all;
        end
        
        % ======================================================================
        %> @brief Return responses for all nodes and an interpolated keff.
        %> @param   keff	Value of keff to interpolate responses.
        %> @return          Responses    
        % ======================================================================    
        function [R, F, A, L] = get_responses(this, keff)

            dim = this.d_dimension;
            len_R = 2*dim * ...
                    this.d_number_groups *  ...
                    (1+this.d_max_order_space) * ...
                    (1+this.d_max_order_azimuth) * ...
                    (1+this.d_max_order_polar);
              
            R_all = cell(this.d_number_nodes, 1);
            F_all = cell(this.d_number_nodes, 1);
            A_all = cell(this.d_number_nodes, 1);
            L_all = cell(this.d_number_nodes, 1);
            % I'm being a bit lazy: if we need low orders, just reduce the full
            % set and interpolate that.  A better way might be to store only the
            % low order k-dependent values, but this depends on purpose.
            if (this.d_order_space   < this.d_max_order_space   || ...
                this.d_order_azimuth < this.d_max_order_azimuth || ...
                this.d_order_polar   < this.d_max_order_polar   )    
                
                for ni = 1:this.d_number_nodes
                    [R_all{ni}, F_all{ni}, A_all{ni}, L_all{ni}] = ...
                        reduce(this, this.d_R_all{ni}, this.d_F_all{ni}, ...
                               this.d_A_all{ni}, this.d_L_all{ni});
                end
                
            else
                R_all = this.d_R_all; 
                F_all = this.d_F_all;
                A_all = this.d_A_all;
                L_all = this.d_L_all;
    
            end
            len_R = length(R_all{1}(:, 1));
            % Interpolated values
            R = zeros(len_R, len_R, this.d_number_nodes);
            F = zeros(len_R,        this.d_number_nodes);
            A = zeros(len_R,        this.d_number_nodes);
            L = zeros(len_R, 2*dim, this.d_number_nodes);
            
            int = this.d_interp_method;
            
            for ni = 1:this.d_number_nodes
                keffv = this.d_keff_vector{ni};
                nk    = length(keffv);
                %
                tmpR = reshape(R_all{ni}(:, :, :), len_R^2, nk);
                R(:, :, ni) = reshape( ...
                    interp1(keffv, tmpR', keff, int, 'extrap'), len_R, len_R);
                %
                tmpF(:, :) = F_all{ni}(:, :)';
                F(:, ni) = interp1(keffv, tmpF, keff, int, 'extrap')';
                %
                tmpA(:, :) = A_all{ni}(:, :)';
                A(:, ni) = interp1(keffv, tmpA, keff, int, 'extrap')';     
                %
                tmpL = reshape(L_all{ni}(:, :, :), len_R*2*dim, nk);
                L(:, :, ni) = reshape( ...
                    interp1(keffv, tmpL', keff, int, 'extrap'), len_R, 2*dim);
            end
            
            % Save this latest evaluation.
            this.d_R = R;
            this.d_F = F;
            this.d_A = A;
            this.d_L = L;

        end        
        
        function set_interp_method(this, method)
            this.d_interp_method = method; 
        end
        
        % ======================================================================
        %> @brief Append a DB to another.
        %>
        %> This assumes that the main header information is the same.
        %> Also, this is for joing HDF5's.  Also, only new assemblies can
        %> be joined, i.e. db's with different orders aren't allowed.
        %>
        %> @param   new_db	ResponseDB object to add to this.
        % ======================================================================            
        function append(this, new_db)
            % Read the new db if not read.
            new_db.read_response();

            new_number_nodes = new_db.number_nodes();
            
            % Append the keff vectors and node descriptions.
            keffs = new_db.keff_vector();
            descs = new_db.node_descriptions();            
            for i = 1: new_number_nodes
                this.d_keff_vector{this.d_number_nodes+i} = keffs{i};
                this.d_node_descriptions{this.d_number_nodes+i} = descs{i};
            end
         
            % Get the responseses.
            [R, F, A, L] = new_db.get_all_responses();
            
            for i = 1:new_number_nodes
                node_index = i + this.d_number_nodes;
                write_response(this, ...
                               node_index, ...
                               this.d_node_descriptions{node_index}, ...
                               R{i}(:,:,:), ...
                               F{i}(:,:), ...
                               A{i}(:,:), ...
                               L{i}(:,:,:));
                this.d_R_all{node_index} = R{i};
                this.d_F_all{node_index} = F{i};
                this.d_A_all{node_index} = A{i};
                this.d_L_all{node_index} = L{i};
            end
            
            % Update number of nodes
            this.d_number_nodes = this.d_number_nodes + new_number_nodes;
            h5writeatt(this.d_file, '/', 'number_nodes', this.d_number_nodes);
            
        end
        
        function k = keff_vector(this)
            k = this.d_keff_vector;
        end
        
        function n = node_descriptions(this)
            n = this.d_node_descriptions;
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
        