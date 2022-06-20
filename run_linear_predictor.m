
function signal_est = run_linear_predictor(signal_with_wgn, samples, order)
it = floor(length(signal_with_wgn)/samples);
signal_est = zeros(it*samples, 1);
for i = 1:it
    S = signal_with_wgn(1+samples*(i-1):samples*i);

    h = linear_predictor(S, order);
    
    signal_est(1+(i-1)*samples:i*samples) = estimate_lp(S, h);
end