# COE352-FEM-Project

(1) Solve first by using a forward Euler time derivative discretization with a time-step of dt = 1/551. Plot the results at the final time. Increase the time-step until you find the instability. What dt does this occur at? How does the solution change as N decreases?

At a time step of dt = 1/551, our results are unstable, as you can see from the plot below. If we increase the timestep, as demonstrated with the plots of dt = 1/500 and 1/400, the results remain unstable and do not converge to the analytical solution. 

![plot][./plots/Forward551.png]
![plot][./plots/Forward500.png]
![plot][./plots/Forward400.png]

However, if we decrease the 
