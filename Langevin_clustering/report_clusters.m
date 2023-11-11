% Write the 1xn sequence of clusters to a single-line file
function[] = report_clusters(clusters, filename, sep)
    % Open a file for writing (replace 'filename.txt' with your desired file name)
    fileID = fopen(filename, 'w');
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Error opening file for writing.');
    end
    
    % Write the vector to the file as a single line, separating each element with a comma
    if isnan(sep)
        fprintf(fileID, '%d', clusters(1:end));
    else
        writematrix(clusters, filename, 'Delimiter', sep);
    end
    %fprintf(fileID, '%d', clusters(1:end));
    %fprintf(fileID, '%d\n', clusters(end));  % New line at the end
    
    % Close the file
    fclose(fileID);
end