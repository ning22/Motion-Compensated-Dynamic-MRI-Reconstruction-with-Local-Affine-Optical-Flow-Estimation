function w = cw2separate(w, J)
% 2nd part of cplxdual2D after icplxdual2D_sym
for j = 1:J+1
    if j == J+1    
        [wp_r wm_r] = pm(w{j}{1}{1},w{j}{1}{2});
        [wp_i wm_i] = pm(w{j}{2}{1},w{j}{2}{2});
        w{j}{1}{1} = wp_r;
        w{j}{2}{2} = -wm_r;
        w{j}{2}{1} = -wp_i;
        w{j}{1}{2} = -wm_i;
    else
        for m = 1:3
            %         [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
            %         [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
            %
            % I think it should have been
            switch m
                case 1
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{1}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{2}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{2}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = -wp_i;
                    w{j}{1}{2}{m} = wm_i;
                case 2
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{1}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{2}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{2}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = wp_i;
                    w{j}{1}{2}{m} = -wm_i;
                case 3
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{1}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{2}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{2}{2}{m} = -wm_r;
                    w{j}{2}{1}{m} = wp_i;
                    w{j}{1}{2}{m} = wm_i;
                    
                otherwise
                    disp('na re na');
            end
        end
    end
end
