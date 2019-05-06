# mitgcm_configs

## Barotropic eddy over fine scale topography

![eddy over topography](input/eddy-iwave.png)

## Time evolution

we set up the MITgcm
to simulate a barotropic eddy (20km radius) over a sinusoidal
topography ($\lambda = 1$ km).

We set the stratification to $N =
10^{-3}$ s$^{-1}$. The eddy is initially in geostrophic balance: the
pressure anomaly comes exclusively from a sea surface height
anomaly. In cylindrical coordinates centered on the eddy axis of
rotation, we have

$$
fu_\theta = -\frac{\partial P}{\partial r}\, ,
$$

with $u_\theta$ the azimuthal velocity and $r$ the radial
coordinate. Here and in the following equations, we assume that the
eddy is axisymmetric.

After one day, the velocity profile of the eddy changes
dramatically with a strong damping in the lower kilometer (just above
the topography), and we were wondering how does this damping
occur. The linear equation for the azimuthal velocity is

$$
\frac{\partial u_\theta}{\partial t} = -f u_r + \nu \frac{\partial^2 u_\theta}{\partial z^2} + \frac{\partial }{\partial z} \overlin\
e{u'_\theta w} \, ,
$$

with the last term $\overline{u'_\theta w}$ the Reynolds stress (where
the overline is for the azimuthal average) and we omitted the pressure
gradient (because of the axisymmetric constraint). Among the three
terms in the rhs of the above equation, the Coriolis term is the main
term responsible for the eddy 'damping'. The two other terms are still
keys because they are at the origin of the process: in fact, because
the Reynolds stress is non zero, it will decrease $u_\theta$ which
will disrupt the geostrophic balance and trigger and inward flow (for
an anticylone) above the topography. Back in the equation of evolution
of $u_\theta$ we see that this inward flow will force $u_\theta$ to
grow. Since $u_\theta$ is negative for an anticyclone, it looks like a
damping of the eddy but is in fact a gesotrophic adjustment (angular
momentum redistribution). In the three figures below, we plotted a
vertical cross section of the eddy (azimuthally averaged). In all
plots, the contour lines are the azimuthal velocity after one day of
the evolution (dashed line: negative velocity for an anticyclone). The
first plot is the time derivative of the azimuthal velocity, the
second plot is the Reynolds stress term, and the last plot is the
coriolis term.

![](dudt.png)
![](dupwpdz.png)
![](ucori.png)

