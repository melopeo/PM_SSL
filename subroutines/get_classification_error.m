function error_vec = get_classification_error(C, GroundTruth, idxSample)

    yu               = GroundTruth;
    yu(idxSample)    = [];
    
    y_hat            = C;
    y_hat(idxSample) = [];
    
    error_vec        = mean( y_hat ~= yu );