function ipdg = build_method_flags(ipdg)
%BUILD_METHOD_FLAGS Set beta according to IPDG method.

    method_name = upper(ipdg.method);

    switch method_name
        case 'SIPG'
            ipdg.beta = 1;
        case 'NIPG'
            ipdg.beta = -1;
        case 'IIPG'
            ipdg.beta = 0;
        otherwise
            error('Unknown IPDG method: %s', ipdg.method);
    end
end