function [dataset, dataset_outlier, dataset_outlier_noise] = create_dataset(n)     
    for D = [4 9 30]
        k = 5;
        N1 = 70 * D;
        for i = 1:k
            b = normc(randn(D, 1));
            U = null(b');
            P = U * U';
            X = P * normc(randn(D, N1));

            X6 = X;
            X8 = X;
            X10 = X;
            for nn = 2:n
                b = normc(randn(D, 1));
                U = null(b');
                P = U * U';
                Y10 = P * normc(randn(D, 1.0^(nn-1)*N1));
                Y6 = Y10(:, 1:ceil(0.6^(nn-1)*N1));
                Y8 = Y10(:, 1:ceil(0.8^(nn-1)*N1));

                eval(['dataset.D', num2str(D), '_n', num2str(nn), '_alp6_', num2str(i), '=[X6 Y6];'])
                eval(['dataset.D', num2str(D), '_n', num2str(nn), '_alp8_', num2str(i), '=[X8 Y8];'])
                eval(['dataset.D', num2str(D), '_n', num2str(nn), '_alp10_', num2str(i), '=[X10 Y10];'])
                X6 = [X6 Y6];
                X8 = [X8 Y8];
                X10 = [X10 Y10];
            end
        end
    end
    
    
%     parms = rand_parms_spec();
%     parms.D = D;
%     parms.n = n;
%     parms.alpha = alp;
%     parms.r = 0;
%     parms.sigma = 0;
%     for i = 1:k
%         info = info_compute(parms);
%         X = info.Xtilde;
%         eval(['dataset.D', num2str(D), '_n', num2str(n), '_alp', num2str(10*alp), '_', num2str(i), '=X;'])
%         M =  ceil(0.1 * info.N / (1 - 0.1));
%         O = normc(randn(parms.D, M));
%         eval(['dataset_outlier.D', num2str(D), '_n', num2str(n), '_alp', num2str(10*alp), '_', num2str(i),'=[X, O];'])
%         X = normc(X + 0.01 * randn(D, size(X, 2)));
%         eval(['dataset_outlier_noise.D', num2str(D), '_n', num2str(n), '_alp', num2str(10*alp), '_', num2str(i),'=[X, O];'])
%     end
            
end

