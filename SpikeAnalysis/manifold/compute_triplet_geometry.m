
        function [proj, perp, adist] = compute_triplet_geometry(R, anc, mid, post)
            n = length(anc);
            proj  = zeros(n,1);
            perp  = zeros(n,1);
            adist = zeros(n,1);
            for k = 1:n
                v_a  = R(anc(k),  :);
                v_m  = R(mid(k),  :);
                v_p  = R(post(k), :);

                ax   = v_m - v_a;
                adist(k) = norm(ax);
                ax   = ax / (norm(ax) + 1e-12);   % unit anchor→middle

                disp = v_p - v_a;                  % post - anchor
                proj(k) = dot(disp, ax);
                perp(k) = norm(disp - proj(k)*ax);
            end
        end