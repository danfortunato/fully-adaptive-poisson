function dom = kdrefine(dom, opts)

arguments
    dom
    opts.fac = 1/3;
end

n = dom.N;
done = false;
while ( ~done )

    kd = KDTree(cell2mat(dom.x));

    np = dom.np;
    split = false(np, 1);
    for k = 1:np
        % Exclude nodes from self and neighboring panels:
        kl = mod((k-1)-1, np) + 1;
        kr = mod((k+1)-1, np) + 1;
        iexclude = [n*(kl-1)+(1:n), n*(k-1)+(1:n), n*(kr-1)+(1:n)].';
        exclude = false(n*np, 1);
        exclude(iexclude) = true;

        % Crudely find near points by doing two kd-ball queries centered at the
        % first and last panel nodes:
        arclen = sum(dom.w{k});
        r = opts.fac * arclen;
        inear1 = kd.ball(dom.x{k}(1,:), r);
        inear2 = kd.ball(dom.x{k}(n,:), r);
        near = false(n*np, 1);
        near(inear1) = true;
        near(inear2) = true;

        % Split this panel if the near field includes any non-neighboring panels
        if ( any(near & ~exclude) )
            split(k) = true;
        end
    end

    if ( any(split) )
        ids = find(split);
        dom = refine(dom, 1, ids);
    else
        done = true;
    end

end

end
