function casted = type_cast(data, type)

if strcmp(type, 'uint8')
    casted = uint8(floor(data*255));
elseif strcmp(type, 'uint16')
    casted = uint16(floor(data*65536));
elseif strcmp(type, 'single')
    casted = single(data);
elseif strcmp(type, 'double')
    casted = data;
end