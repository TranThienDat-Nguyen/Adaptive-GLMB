function state = meas2state(z)
state = zeros(4,1) ; 
state([1,3],1) = z ; 