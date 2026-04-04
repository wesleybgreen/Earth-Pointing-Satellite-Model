# Earth pointing satellite model

This is a personal project I've been working on in my spare time, while pursuing my MS with thesis at Virginia Tech. This is an attitude pointing model that would take place after a detumbling mode, where the model stays pointed towards the Earth for an amount of time. I utilize a multiplicative EKF with a sun sensor, magnetometer, and gyroscopes for attitude control. I also technically have a star tracker that outputs a quaternion directly; however, I decided to model these other sensors to better test the algorithm and because sometimes a star tracker isn't always used. The only issues that arise with the attitude control are during eclipse, when it is a little bit shaky; however, I may look into an Earth-horizon sensor instead of the sun sensor to deal with that. The output plot for the attitude error can be seen below, where it converges at 0 for the quaternion vector and 1 for the scalar, with the MEKF tracking the true values. 

![Quaternion Error](QuatError.jpeg)

The eclipse can be seen to occur during the 3000-4000 second area, where it has a slight hiccup. Next, the angular velocity of the satellite can also be seen to approach 0 below:

![Angular Velocity Error](AngVel.jpeg)

Finally, the gyroscope bias plots can be seen below, where I have zoomed in on the MEKF estimate to show the accuracy compared to the true value and then also the zoomed-out version.

![Gyro Bias Zoom](GyroBiasZoom.jpeg)
![Gyro Bias](GyroBias.jpeg)

The MEKF estimate follows pretty closely to the true value, but not exactly. I plan to look more into this and determine more specific hardware specs for sensor accuracy, rather than just going with benchmark values. I may also experiment with different combos of sensor models just to learn about each kind of sensor, but also, sensor combos are specific to mission needs, so I don't have a specific goal in mind for model accuracy. 

Plans moving forward:

I plan to add a detumbling mode next, which then moves into the pointing mode. After that, I would like to add a momentum desaturation mode since the reaction wheels eventually reach saturation.
