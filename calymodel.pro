function calymodel,x,param
N  = param[0]
fc = param[1]
tau  = x*N/3.768
ycal = 1-fc*(1.-exp(-1.*tau))
return, ycal
end
