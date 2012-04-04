% Include the MTP source directory
path(path, '../')

% Open a DB
input1 = Input();
put(input1, 'rf_db_name',        'c5g7_assemblies.h5');
db1 = ResponseDB(input1);
read_response(db1);

input2 = Input();
put(input2, 'rf_db_name',        'c5g7_assemblies2.h5');
db2 = ResponseDB(input2);
read_response(db2);

input3 = Input();
put(input3, 'rf_db_name',        'c5g7_assemblies3.h5');
db3 = ResponseDB(input3);
read_response(db3);

db1.append(db2);
db1.append(db3);
