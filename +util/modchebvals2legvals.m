function legvals = modchebvals2legvals(chebvals)

chebvals = util.modchebvals2chebvals(chebvals);
legvals = chebvals2legvals(chebvals);

end
