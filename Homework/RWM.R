## Random Walk Metropolis ##

# inputs:
#   q_cur: initial state of q, lo<q<up
#   u_cur: initial potential energy
#   U: =-log(density(q)), potential function of q
#   eps: step size
# outputs:
#   q: new state
#   u: new potential energy
#   Ind: proposal acceptance indicator

RWM = function (q_cur, u_cur, U, eps=.1){
#    browser()
    # initialization
    q = q_cur; D = length(q)
	u = u_cur;
    
    # sample velocity
    z = rnorm(D) # standard multi-normal
	
    # evolve position using random direction
	q = q + sqrt(eps)*z
	
	# Evaluate new potential energy
	u = U(q)
    
    # Accept or reject the state and check constraint, returning either
    # the position at the end of the trajectory or the initial position
    logAP = -u + u_cur
	
	if( is.finite(logAP)&&(log(runif(1))<min(0,logAP)) ){
		return (list(q = q, u = u, Ind = 1))
	}else return (list(q = q_cur, u = u_cur, Ind = 0))
    
}
