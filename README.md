# Earth pointing satellite model

This is a personal project I have been working on that would take place post-detumbling. It simulates the pointing of a satellite
that is pointed towards the Earth over time. It involves a "perfect sim", which is what I started with and just covers the dynamics.
From there, I have since added a multiplicative EKF and several sensors (involving synthetic sensor measurements), with plans to
add more modes, such as detumbling, soon. To run this, you need to have numpy, scipy.linalg, scipy.integrate, matplotlib.pyplot, and
tqdm installed. 
