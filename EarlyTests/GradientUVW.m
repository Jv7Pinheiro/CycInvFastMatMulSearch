function GradientUVW(term)
    % Open the file for reading
    fileID = fopen('ListUVW_Output.txt', 'r'); % Replace 'file_name.txt' with the file path/name

    if fileID == -1
        error('File not found or could not be opened.');
    end

    % Read the whole file into a string
    fullText = fscanf(fileID, '%c');
    fclose(fileID);

    % Split the text into lines
    lines = strsplit(fullText, '\n');

    % Process each line separately to find the term
    for lineNum = 1:numel(lines)
        line = lines{lineNum};

        % Find all occurrences of the term in the line
        indices = strfind(line, term);

        % Process each occurrence and find the "word"
        for i = 1:numel(indices)
            startIndex = indices(i);

            % Find the previous space character
            prevSpace = find(line(1:startIndex-1) == ' ', 1, 'last');

            % Find the next space character
            nextSpace = find(line(startIndex:end) == ' ', 1, 'first') + startIndex - 1;

            % Extract the "word" that contains the term
            word = line(prevSpace+1:nextSpace-1);

            % Output the line reference and word to the Command Window
            %fprintf('(Line %d) : %s\n', lineNum, word);
            fprintf('%s : %s\n', line(1:3), word);
        end
    end
end