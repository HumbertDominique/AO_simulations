%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 01.07.2026 - Dominique Humbert - Submition version
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cfg = readConfig(filename)
% READCONFIG  Read a key=value parameter file into a struct.
%
%   cfg = readConfig('ao_inputs.txt')
%
%   Lines beginning with '#' and blank lines are ignored.


    cfg = struct();
    fid = fopen(filename, 'r');
    if fid == -1
        error('readConfig:fileNotFound', 'Cannot open config file: %s', filename);
    end

    while ~feof(fid)
        raw = fgetl(fid);
        if ~ischar(raw), break; end

        % Strip inline comments and whitespace
        commentIdx = strfind(raw, '#');
        if ~isempty(commentIdx)
            raw = raw(1:commentIdx(1)-1);
        end
        line = strtrim(raw);
        if isempty(line), continue; end

        eqIdx = strfind(line, '=');
        if isempty(eqIdx), continue; end

        key = strtrim(line(1:eqIdx(1)-1));
        val = strtrim(line(eqIdx(1)+1:end));
        if isempty(key) || isempty(val), continue; end

        % Parse value
        if ~isempty(val) && val(1) == '[' && val(end) == ']'
            cfg.(key) = str2num(val);
        elseif strcmpi(val, 'true')
            cfg.(key) = true;
        elseif strcmpi(val, 'false')
            cfg.(key) = false;
        else
            num = str2double(val);
            if ~isnan(num)
                cfg.(key) = num;
            else
                cfg.(key) = val;  
            end
        end
    end

    fclose(fid);
end
