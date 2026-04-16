function hash = DataHash(data)
% DATAHASH - Reliable MD5 hash for SVD cache invalidation.
% Uses Java MessageDigest (available in all MATLAB versions with JVM).
% For complex data, hashes the interleaved real+imag byte stream.

if isempty(data)
    hash = 'empty';
    return;
end

try
    import java.security.MessageDigest;
    md = MessageDigest.getInstance('MD5');

    if isnumeric(data) || islogical(data)
        d = double(data(:));
        if ~isreal(d)
            d = [real(d); imag(d)];
        end
        bytes = typecast(d, 'uint8');
    elseif ischar(data) || isstring(data)
        bytes = uint8(char(data(:))');
    else
        bytes = uint8(sprintf('%s_%s', class(data), mat2str(size(data))));
    end

    md.update(bytes);
    digest = md.digest();
    hash   = sprintf('%02x', typecast(digest, 'uint8'));

catch
    % Fallback if JVM unavailable (e.g., MATLAB compiled mode)
    sz = size(data);
    if isnumeric(data) && numel(data) > 0
        hash = sprintf('%s_%s_%.10g_%.10g_%.10g', mat2str(sz), class(data), ...
            double(data(1)), double(data(end)), double(sum(data(:))));
    else
        hash = sprintf('%s_%s', class(data), mat2str(sz));
    end
end

end