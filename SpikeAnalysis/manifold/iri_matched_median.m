
function [mu_st, mu_sp, perm_p] = iri_matched_median(metric, ...
                                      kp_st, kp_sp, dt_vec, nPerm)

    dt_bins  = [0 logspace(log10(0.05), log10(30), 8) Inf];
    st_idx   = find(kp_st);
    sp_idx   = find(kp_sp);
    dt_bin   = discretize(abs(dt_vec), dt_bins);

    matched_st = [];
    matched_sp = [];

    for ii = 1:length(st_idx)
        b     = dt_bin(st_idx(ii));
        cands = sp_idx(dt_bin(sp_idx) == b);
        if isempty(cands), continue; end
        pick  = cands(randi(length(cands)));
        matched_st = [matched_st; metric(st_idx(ii))];
        matched_sp = [matched_sp; metric(pick)];
    end

    if length(matched_st) < 5
        mu_st = nan; mu_sp = nan; perm_p = nan; return
    end

    mu_st = nanmean(matched_st);
    mu_sp = nanmean(matched_sp);
    obs   = mu_st - mu_sp;

    pooled = [matched_st; matched_sp];
    n      = length(matched_st);
    perms  = zeros(nPerm,1);
    for pp = 1:nPerm
        idx      = randperm(length(pooled));
        perms(pp) = nanmean(pooled(idx(1:n))) - ...
                    nanmean(pooled(idx(n+1:end)));
    end

    perm_p = mean(abs(perms) >= abs(obs));

end