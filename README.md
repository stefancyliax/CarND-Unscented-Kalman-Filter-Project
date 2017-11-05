# Udacity Project 7: Unscented Kalman Filter
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

The goal of this project was to build on the first project but now use a more sophisticated Kalman Filter. The Unscented Kalman Filter. The project uses the same LIDAR and RADAR data as the first one.

In the end I was able to achieve a RMSE of:

| State | RSME of Unscented Kalman Filter | RMSE of Extended Kalman Filter |
|-------|---------------------------------|--------------------------------|
| X     | __0.0701__                          | 0.0973                         |
| Y     | __0.0846__                          | 0.0855                         |
| VX    | __0.1710__                          | 0.4513                         |
| VY    | __0.2007__                          | 0.4399                         |

Project video: (Youtube link)

[![Project track](https://github.com/stefancyliax/CarND-Unscented-Kalman-Filter-Project/raw/master/output_images/project_video.gif)](https://youtu.be/eaKn6NKWmrA)

### Basic principle

The basic working principle of a Kalman filter is to track the state of an object using state and the uncertainty of this state. Then, at some point, a new sensor measurement arrives, again with a certain uncertainty.
To update our state with the new measurement, we first do a prediction of what the state should be like at the given point in time. Then we weight the prediction and the sensor measurement by their distinct uncertainty. If the measurement is very uncertain, the predicted state is weighted higher and vise versa.
This allows for an algorithm that combines several noisy measurements and arrive gradually at an accurate state.

### Project Approach

For P2 only the structure of the program were given. Because of the complex process of the Unscented Kalman Filter, quite a bit of programming and debugging were involved. P2 used the same [simulator](https://github.com/udacity/self-driving-car-sim/releases) as P1.


### Basic Build Instructions

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

For dependencies and setup see the Udacity starter [documentation](https://github.com/stefancyliax/CarND-Unscented-Kalman-Filter-Project/blob/master/Starter_code_README.md). 
