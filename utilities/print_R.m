fileID = fopen('block_big.txt', 'w');

for i = 1:length(R)
    for j = 1:length(R)
        fprintf(fileID,'%12.8f\n', R(i, j));
    end
end
