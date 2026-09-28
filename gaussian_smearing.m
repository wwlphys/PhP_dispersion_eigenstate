function val=gaussian_smearing(x,x0,sigma)
  val=1.0/(sigma*sqrt(2*pi))*exp(-((x-x0)^2)/(2*(sigma^2)));
end 
