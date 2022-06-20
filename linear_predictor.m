%Devuelve los coeficients del predictor según el orden y el vector de muestras.
function h = linear_predictor(S, order)
    [normalizedACF, lags] = autocorr(S, numLags=order);
    unnormalizedACF = normalizedACF * var(S);
    
    
    R_SS = zeros(order);
    %r_SS = zeros(order);
    
    for i = 1:order
        for j = 1:order
            R_SS(i,j) = unnormalizedACF(abs(i-j)+1);
            %r_SS(i,j) = normalizedACF(abs(i-j)+1);
        end
    end
    
    R_SS_vector = zeros(order, 1);
    %r_SS_vector = zeros(order, 1);
    for i = 1:order
        R_SS_vector(i) = unnormalizedACF(i+1);
        %r_SS_vector(i) = normalizedACF(i+1);
    end
    
    h = R_SS \ R_SS_vector;
    %h2 = r_SS \ r_SS_vector; % h = h2
end