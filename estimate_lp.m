% Estimo S
function Sest = estimate_lp(Swgn, h)
    samples = length(Swgn);
    order = length(h);
    
    Sest = zeros(samples,1);
    Sest(1:order) = Swgn(1:order);
    for k = order+1:samples
        Sest(k) = h' * Swgn(k-order:k-1);
    end
end