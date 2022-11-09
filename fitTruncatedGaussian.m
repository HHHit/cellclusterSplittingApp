function [varall] = fitTruncatedGaussian(input, ratio)
    input = input - mean(input);
    signinput = sign(input);
    [y,sid] = sort(abs(input));
    signinput = signinput(sid);
    l = length(input);
    if l < 20
        varall = var(input);
    else
        rmc = ceil(l*ratio);
        out = y(1: (l - rmc+1)).*signinput(1: (l - rmc+1));
        ratio_left = sum(signinput(l - rmc:l) == -1)/l;
        ratio_right = sum(signinput(l - rmc:l) == 1)/l;
    if(ratio_left ~= 0 || ratio_right~=0)
        if(ratio_left == 0)
            beta = norminv(1 - ratio_right);
            phi_beta = normpdf(beta);
            var_ratio = 1 + (- beta*phi_beta)/(1 - ratio_right) ...
            - ((phi_beta)/(1 - ratio_right))^2;
            varall = var(out)/var_ratio;
        elseif(ratio_right == 0)
            alpha = norminv(ratio_left);
            phi_alpha =  normpdf(alpha);
            var_ratio = 1 + (alpha*phi_alpha)/(1 - ratio_left) ...
            - ((phi_alpha)/(1 - ratio_left))^2;
            varall = var(out)/var_ratio;
        else
        alpha = norminv(ratio_left);
        beta = norminv(1 - ratio_right);
        phi_alpha =  normpdf(alpha);
        phi_beta = normpdf(beta);
        var_ratio = 1 + (alpha*phi_alpha - beta*phi_beta)/(1 - ratio_left - ratio_right) ...
        - ((phi_alpha - phi_beta)/(1 - ratio_left - ratio_right))^2;
        varall = var(out)/var_ratio;
        end
    else
        varall = var(input);
    end
    end
end
