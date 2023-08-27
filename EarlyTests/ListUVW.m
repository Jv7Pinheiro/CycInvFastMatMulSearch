function ListUVW(length_ijk, length_l)
    % Replace 'file_name.txt' with the desired file name and path.
    fid = fopen('ListUVW_Output.txt', 'w'); % 'w' opens the file in write mode.
    
    fprintf(fid, "(i, j, k) -> [%d], (l) -> [%d]\n\n\n", length_ijk, length_l);
    for i = 1:length_ijk
        for j = 1:length_ijk
            for k = 1:length_ijk
                fprintf(fid, "%d%d%d : ", i, j, k);
                for l = 1:length_l
                    fprintf(fid, "( U%d%d*V%d%d*W%d%d ", i, l, j, l, k, l);
                    fprintf(fid, "+ W%d%d*U%d%d*V%d%d ", i, l, j, l, k, l);
                    fprintf(fid, "+ V%d%d*W%d%d*U%d%d ) + ", i, l, j, l, k, l);
                end
                fprintf(fid, "\n");
            end
        end
    end
    fprintf("Done\n")
    fclose(fid);
end
