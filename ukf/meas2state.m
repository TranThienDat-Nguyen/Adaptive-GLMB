function state = meas2state(z)
state = zeros(5,1) ; 
state(1) = z(2) * sin(z(1)) ; 
state(3) = z(2) * cos(z(1)) ; 