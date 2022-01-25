function X = scale_interval(X, dom1, dom2)

X = (X - dom1(1)) ./ diff(dom1); % Map dom1 to [0 1]
X = X*diff(dom2) + dom2(1);      % Map [0 1] to dom2

end
