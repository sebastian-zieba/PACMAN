.. _models:

fit_par.txt and Fitting models
================================

Let's have a look at the currently implemented models in pacman:

Models
--------------------------------

Instrument Systematics
''''''''''''''''''''''''''''''

* `constant.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/constant.html#constant>`_

  free parameters: c

.. note:: c is in log10. A average flux of 10^7 therefore leads to approximately c = 7.

* `model_ramp.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/model_ramp.html#model_ramp>`_

  free parameters: r1, r2, r3

* `polynomial1.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/polynomial1.html#polynomial1>`_

  free parameters: v

* `polynomial2.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/polynomial2.html#polynomial2>`_

  free parameters: v, v2

* `upstream_downstream.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/upstream_downstream.html#upstream_downstream>`_

  free parameters: scale

* `divide_white.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/divide_white.html#divide_white>`_



Astrophysical
''''''''''''''''''''''''''''''

* `transit.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/transit.html#transit>`_

  free parameters: t0, per, rp, a, inc, ecc, w, u1, u2, limb_dark

* `eclipse.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/eclipse.html#eclipse>`_

  free parameters: t_secondary, per, rp, fp, a, inc, ecc, w

* `sine1.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine1.html#sine1>`_

  free parameters: a1, omega1, phi1

* `sine2.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine2.html#sine2>`_

  free parameters: a1, omega1, phi1, a2, omega2, phi2, a3, omega3, phi3, a12, omega12, phi12, a22, omega22, phi22, a32, omega32, phi32

* `sine_curve.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine_curve.html#sine_curve>`_

  free parameters: amp1, theta1, per, amp2, theta2


fit_par.txt file
----------------------------

This file has to be set up when running Stage 30.
Here's an example:

.. include:: media/fit_par.txt
   :literal:

Let's have a look at each column:

- parameter

  A list of the different parameters in each model is above.


- fixed

  If set to False, the parameter will be a free parameter.


- tied

  If the user wants to tie a parameter over all visits, set -1.

  If the user does not want to tie a certain parameter, he or she has to duplicate the line as often as they have visits.

  Example: c in the template above. The code assumes that the user sorted the rows in the correct order.


- value

  If fixed was set to True, this will be the used value for the parameter.

  If fixed was set to False, this is the initial guess for the parameter.


- lo_lim

  Use lower bounds for the least squares routine?


- lo_val

  lower bound value for the least squares routine.


- hi_lim

  Use upper bounds for the least squares routine?


- hi_val

  upper bound value for the least squares routine.


- prior

  Prior for the sampling?

  * X: No prior

  * U: uniform prior

  * N: Gaussian prior


- p1 & p2

  If prior = U -> lower and upper bounds for the uniform prior

  if prior = N -> mean and 1 sigma for the gaussian prior


- step_size

  Sets a step_size for the least squares and sampling.
