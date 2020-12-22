# quadrotor-state-estimator
 
Simulation of quadrotor state estimation using **Complementary Filter**, **Extended Kalman Filter** and **Error-state Kalman Filter** written in MATLAB.

<img src="https://github.com/shengwen-tw/quadrotor-state-estimator/blob/master/images/3d-trajectory.png?raw=true" width="40%" height="40%">

## Dataset

The CSV file is consist of the sensor measurements from **accelerometer**, **gyroscope**, **magnetometer**, **gps receiver** and **barometer**

## Run simulation

**Octave under GNU/Linux**

```
> octave main.m
```

**MATLAB**

```
execute main.m
```

## References

[1] [Quaternion kinematics for the error-state Kalman flter](http://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf)

[2] [A Double-Stage Kalman Filter for Orientation Tracking With an Integrated Processor in 9-D IMU](https://ieeexplore.ieee.org/document/6316172)

[3] [Keeping a Good Attitude: A Quaternion-Based Orientation Filter for IMUs and MARGs](https://www.mdpi.com/1424-8220/15/8/19302/pdf)
