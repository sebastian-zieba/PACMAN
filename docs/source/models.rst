.. _models:

fit_par.txt and Fitting models
================================

Instrument Systematics
''''''''''''''''''''''''''''''

* `constant.py <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/lib/models/constant.html#constant>`_

free parameters: c

.. note:: c is in log10. A average flux of 10^7 therefore leads to c = 7.

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
