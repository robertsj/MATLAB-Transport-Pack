% Include the MTP source directory
path(path, '../')

% Open a DB from hdf5 file
input1 = Input();
put(input1, 'rf_db_name',        'c5g7_assemblies.h5');
put(input1, 'rf_db_name_mat',    'c5g7_assemblies.mat');

db1 = ResponseDB(input1);
read_response(db1);

% Write it back out as mat file
[R, F, A, L] = db1.get_all_responses();
nd =  db1.node_descriptions();
db1.initialize_write_mat();
for node_index = 1:length(nd)
    db1.write_response_mat(...
        node_index, ...
        nd{node_index}, ...
        R{node_index}(:,:,:), ...
        F{node_index}(:,:), ...
        A{node_index}(:,:), ...
        L{node_index}(:,:,:));
end

% db2 = ResponseDB(input1);
% read_response_mat(db2);
% [R2, F2, A2, L2] = db2.get_all_responses();