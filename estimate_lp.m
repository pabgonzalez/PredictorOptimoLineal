% Estimo S
function Sest = estimate_lp(Swgn, h)
    samples = length(Swgn);
    order = length(h);
    
    %Agrego ceros al principio para estimar la señal
    Swgn_prima = [zeros(order,1); Swgn];
    Sest = zeros(samples,1);
    for k = 1:samples
        Sest(k) = h' * Swgn_prima(k:k+order-1);
    end
end