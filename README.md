# Simulated Annealing
This is an implementation of the simulated annealing algorithm for Bocconi's 20602 COMPUTER SCIENCE (ALGORITHMS) course.

# Data
There are two data generating processes for this COP. First, there is `data_generator_cirlce` which generators data points on an ellipse.
This can be used for model validation as the shortest path of an ellipse is known. 
Secondly, there is the `random_data_generator` which randomly generates points in a grid. 

# Arguments
The algorithm takes the following arguments:
- **x0**: data points
- **step_max**: maximum step size of the algorithm
- **t_min**: minimum temperatur that can be reached
- **t_max**: initial max temperatur
- **cooling_schedule_type**: Cooling schedule. Implemented are linear, geometric, linear multiplicative and quadratic multiplicative 
- **alpha**: Control parameter in cooling schedules
- **tau**: Number of MCMC runs at each temperature
