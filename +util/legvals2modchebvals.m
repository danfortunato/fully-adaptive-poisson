function chebvals = legvals2modchebvals(legvals)

chebvals = legvals2chebvals(legvals);
chebvals = util.chebvals2modchebvals(chebvals);

end
