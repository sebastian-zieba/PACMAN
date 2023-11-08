.. _models:

Fitting the Light Curves
================================
``PACMAN`` can fit time series observations of exoplanet transits, eclipses, and phase curves. It includes models for both the astrophysical signal and instrument systematic noise (e.g. ramps and slopes).

Here we describe the currently implemented models. Initial guesses and priors for the model parameters are listed in the `fit_par.txt file <https://pacmandocs.readthedocs.io/en/latest/models.html#fit-par-txt-file>`_.

The full model, which will be fit to the dataset, is the product of the physical and systematic flux. This ensures that we fit for the physical and systematic components simultaneously. For an example see `Kreidberg et al., 2018 <https://ui.adsabs.harvard.edu/abs/2018AJ....156...17K/abstract>`_.

Models
--------------------------------

Instrument Systematics
''''''''''''''''''''''''''''''

* `constant.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/constant.html#constant>`_

  A constant (free parameters: c)

.. note:: c is in log10. An average flux of 10^7 photoelectrons therefore leads to approximately c = 7.

* `model_ramp.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/model_ramp.html#model_ramp>`_

  Exponential ramps fit to each orbit (free parameters: r1, r2, r3)

* `polynomial1.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/polynomial1.html#polynomial1>`_

  Linear slope (free parameters: v)

* `polynomial2.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/polynomial2.html#polynomial2>`_

  Quadratic trend (free parameters: v, v2)

* `exponential_visit.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/exponential_visit.html#exponential_visit>`_

  Exponential trend over the whole visit (free parameters: exp1, exp12)

* `logarithmic_visit.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/logarithmic_visit.html#logarithmic_visit>`_

  Logarithmic trend over the whole visit (free parameters: log1, log2)

* `upstream_downstream.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/upstream_downstream.html#upstream_downstream>`_

  Offset between scan directions due to the upstream-downstream effect (free parameters: scale)

* `divide_white.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/divide_white.html#divide_white>`_

  Uses the divide-white method, which assumes that the systematic parameters for the spectroscopic light curves are the same (have the same shape) as for the white light curve. See equation 2 in `Kreidberg et al. (2014) <https://arxiv.org/pdf/1401.0022.pdf>`_ for a reference. This model does not have any additional free parameters (so nothing has to be added to the fit_par.txt file that used).

* `constants_cj.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/constants_cj.html#constants_cj>`_

  Alternative for model_ramp which fits a constant to every j-th exposure in an orbit. See equation 1 in `Kreidberg et al. (2019) <https://arxiv.org/pdf/1904.10618.pdf>`_ and references within for an application of this model. 

.. note:: c is in log10. An average flux of 10^7 photoelectrons therefore leads to approximately c = 7.


* `uncmulti.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/uncmulti.html#uncmulti>`_

  Scales the errorbars at every iteration of the sampler (free parameter: uncmulti_val). Does not work for the least squares fit only for the samplers, mcmc and dynesty.

.. warning:: With the current version of PACMAN, uncmulti_val has to be entered into fit_par.txt as the last parameter!!
 
  This function has been implemented as an alternative to the 'rescale_uncert' technique. This rescales the errorbars of the flux measurements after the least squares routine so that chi2_red = 1. This might be problematic however, if the least squares is having troubles finding a good solution. Then the errorbars would be overestimated.



Astrophysical
''''''''''''''''''''''''''''''

* `transit.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/transit.html#transit>`_

  Planetary transit (free parameters: t0, per, rp, a, inc, ecc, w, u1, u2).
  The parameters follow the parametrization from batman.
 - t0: transit midtime in days. Note that toffset from the pcf will be subtracted from this parameter to avoid fitting for high numbers.
 - per: period in days
 - rp: planet to star ratio Rp/Rs. unitless.
 - a: semi-major-axis to star radius ratio a/Rs. unitless.
 - inc: inclination in degree
 - ecc: eccentricity. unitless
 - w: argument of periastron in degrees. typically set to 90° if ecc = 0.
 - u1, u2: limb darkening

* `eclipse.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/eclipse.html#eclipse>`_

  Secondary eclipse (free parameters: t_secondary, per, rp, fp, a, inc, ecc, w)
  The parameters follow the parametrization from batman.

 - t_secondary: eclipse midtime in days.
 - per: period in days
 - rp: planet to star ratio Rp/Rs. unitless.
 - fp: planet to star flux Fp/Fs. unitless
 - a: semi-major-axis to star radius ratio a/Rs. unitless.
 - inc: inclination in degree
 - ecc: eccentricity. unitless
 - w: argument of periastron in degrees. typically set to 90° if ecc = 0.

* `sine1.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine1.html#sine1>`_

  Sinusoid (free parameters: a1, omega1, phi1)

* `sine2.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine2.html#sine2>`_

  Sum of three sinusoids (free parameters: a1, omega1, phi1, a2, omega2, phi2, a3, omega3, phi3, a12, omega12, phi12, a22, omega22, phi22, a32, omega32, phi32)

* `sine_curve.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/sine_curve.html#sine_curve>`_

  Sum of two sinusoids (free parameters: amp1, theta1, per, amp2, theta2)


The fit_par.txt file
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



Add a new model
----------------------------

What do I need?
''''''''''''''''''''''''''''''

- A model name

- The parameters of your fitting model


Which files do I have to change?
''''''''''''''''''''''''''''''''''

- src/pacman/lib/formatter.py

- src/pacman/lib/functions.py

- src/pacman/lib/models/MODEL_NAME.py

- src/pacman/data/fit_par.py

What do I have to do?
''''''''''''''''''''''''''''''''''

As an example, we will create a new model which will fit the product of a 
polynomial of 1. order and the upstream downstream effect.
We will call this model ``polynomial1_full``. 

 1) We create a new py file in src/pacman/lib/models/ called polynomial1_full.py


 2) We create the code for this model which combines a fitting for the scale and 
a linear polynomial over the visit:


.. include:: media/models/polynomial1_full.py
   :literal:
 
Note that we use the same function name (``def polynomial1_full(...)``) as file name.


 3) in work.......








