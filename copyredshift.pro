pro copyredshift,science,redshiftdomain

domain = mrdfits('/scr2/nichal/workspace3/SCIENCE/'+redshiftdomain,1)
science.method = redshiftdomain
science.zspec  = domain.zspec
science.zsys_lya = domain.zsys_lya
science.zsys_ism_ave = domain.zsys_ism_ave
science.zsys_ism_std = domain.zsys_ism_std

end
