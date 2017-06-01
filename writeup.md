## Model
My implementation is based on the kinematic model presented in lecture. Below is the
description on each of the component in the model.

### State
My state is a 6-D vector of (x-position, y-position, orientation (psi), velocity along the
orientation (v), cross track error against the y-position (cte), and the orientation error (epsi)).

### Actuators
My actuators contain of two values, which is the steering value and the throttle value.

### Update equations
The update equations for the first 4 states is based on the kinematic vehicle model found in
the lecture, while the update equations for the cross track error and orientation error
is also presented in the lecture.

Refer to line 118 of [MPC.cpp](https://github.com/chunhowt/CarND-MPC-Project/blob/master/src/MPC.cpp#L118)
for more details.

### Cost function
The cost function that the solver tries to optimize consists of three components
- mean square error compared to ideal trajectory (position (cte), orientation (epsi) and
the ideal velocity of the car.
- L2 regularization on the actuator values to prefer minimal usage of steering and accelaration.
- Minimal changes of actuator values between two timestep, to allow smoother driving.

Refer to line 53 of [MPC.cpp](https://github.com/chunhowt/CarND-MPC-Project/blob/master/src/MPC.cpp#L53)
for more details.

### Timestep length and elapsed duration
Here, I picked dt = 0.05 (50 ms) and N = 16 so that we can get close to 1 s of trajectory
to allow the solver to look ahead long enough to find the best trajectory. However,
we didn't make the timestep to be too large since it will be more computationally
expensive to solve for the optimal trajectory.

Initially, I picked dt = 0.1 (100 ms) and N around 12, but noticed that the car will have
problem when we introduce latency to the actuator on the uneven surface (bridge). To get
more fine-grained predictions using similar computation, I halve the dt to 0.05.

## Polynomial fitting and MPC preprocessing
I convert all the waypoints into vehicle coordinates (with orientation equal to
where the vehicle is heading to) before doing the MPC.
To do this, I essentially translate the waypoint location to the vehicle coordinate, and
then do a clockwise rotation on the coordinates. Refer to line 108 of [main.cpp](
https://github.com/chunhowt/CarND-MPC-Project/blob/master/src/main.cpp#L108).

By doing this transformation, it makes the state vector into the MPC much simpler, since
the x,y coordinate are both 0 and the orientation is also zero.

I also need to transform the cte and epsi to be at the vehicle coordinate, see line 118
and 122 of [main.cpp](https://github.com/chunhowt/CarND-MPC-Project/blob/master/src/main.cpp)

## MPC with latency
Initially, I tweak the params by setting latency to 0 ms to minimize the latency effect.
Once I got a working solution with 0 ms actuator latency, I turned it on to 100 ms and realize
that my car will hit the curve on uneven surfaces (eg. bridge).

To solve for this issue, instead of sending the actuator values at timestep 0 to the car,
I sent the actuator value predicted 100ms down the timestep. This helps my car to drive
within the track across the whole track.
