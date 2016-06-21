function tdist, x

nu=8

return, gamma((nu+1)/2)/(sqrt(nu*!pi)*gamma(nu/2))*(1+(x^2)/nu)^(-(nu+1)/2)

end