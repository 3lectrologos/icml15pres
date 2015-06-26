function L = lip(x, y)
  xd = pdist(reshape(x, [], 1));
  yd = pdist(reshape(y, [], 1));
  L = max(yd ./ xd);
end