import warnings
import numpy as np


def _as_bool(value):
    return str(value).lower() == "true"


def prior_bounds(prior, p1, p2, nsigma=5.0):
    """Return lower/upper bounds implied by a prior.

    U: p1, p2 are lower/upper bounds.
    N: p1 is mean, p2 is sigma, bounds are mean +/- nsigma*sigma.
    X: no bounds.
    """
    prior = str(prior).upper()
    p1 = float(p1)
    p2 = float(p2)

    if prior == "U":
        if p2 <= p1:
            raise ValueError(
                f"Uniform prior requires p2 > p1, but got p1={p1}, p2={p2}."
            )
        return p1, p2

    if prior == "N":
        if p2 <= 0:
            raise ValueError(
                f"Normal prior requires sigma p2 > 0, but got p2={p2}."
            )
        return p1 - nsigma * p2, p1 + nsigma * p2

    if prior == "X":
        return -np.inf, np.inf

    raise ValueError(f"Unknown prior type '{prior}'. Expected U, N, or X.")


def validate_fit_par(
    fit_par,
    run_mcmc=False,
    run_nested=False,
    nsigma_lsq=5.0,
    warn_only_for_x=True,
):
    """Validate fit_par.txt.

    - Checks that initial values are inside the implied prior/bounds.
    - Warns if MCMC/nested is requested but a free parameter has prior X.
    """
    problems = []
    x_free = []

    for i in range(len(fit_par)):
        name = str(fit_par["parameter"][i])
        tied = str(fit_par["tied"][i])
        fixed = _as_bool(fit_par["fixed"][i])
        prior = str(fit_par["prior"][i]).upper()
        value = float(fit_par["value"][i])
        p1 = float(fit_par["p1"][i])
        p2 = float(fit_par["p2"][i])

        label = name if tied == "-1" else f"{name}_{tied}"

        if prior not in ["U", "N", "X"]:
            problems.append(
                f"{label}: unknown prior '{prior}'. Use U, N, or X."
            )
            continue

        if (run_mcmc or run_nested) and (not fixed) and prior == "X":
            x_free.append(label)

        if prior in ["U", "N"]:
            lo, hi = prior_bounds(prior, p1, p2, nsigma=nsigma_lsq)

            if not (lo <= value <= hi):
                if prior == "U":
                    problems.append(
                        f"{label}: initial value {value} is outside "
                        f"uniform prior [{lo}, {hi}]."
                    )
                elif prior == "N":
                    problems.append(
                        f"{label}: initial value {value} is outside "
                        f"{nsigma_lsq:g}σ normal-prior sanity range "
                        f"[{lo}, {hi}]."
                    )

    if x_free:
        message = (
            "The following free parameters have prior X while MCMC/nested "
            "sampling is enabled:\n"
            + "\n".join(f"  - {name}" for name in x_free)
            + "\nThese parameters are unconstrained by the sampler prior. "
            "Use prior U for bounded parameters or N for Gaussian priors."
        )

        if warn_only_for_x:
            warnings.warn(message, RuntimeWarning)
        else:
            problems.append(message)

    if problems:
        raise ValueError(
            "\nInvalid fit_par.txt configuration:\n"
            + "\n".join(f"  - {problem}" for problem in problems)
        )


def read_fit_par_for_ls(parinfo, params_s, data, fit_par):
    """Read fit_par.txt and prepare MPFIT parinfo.

    New fit_par format:

    #parameter fixed tied value prior p1 p2

    Bounds:
    - U: [p1, p2]
    - N: [p1 - 5*p2, p1 + 5*p2]
    - X: unbounded; should normally only be used for fixed parameters
    """
    nvisit = data.nvisit
    ii = 0

    for i in range(int(len(data.parnames))):
        tied = str(fit_par["tied"][ii])

        if tied == "-1":
            # One row, shared/tied across visits.
            prior = fit_par["prior"][ii]
            p1 = fit_par["p1"][ii]
            p2 = fit_par["p2"][ii]
            lo, hi = prior_bounds(prior, p1, p2, nsigma=5.0)

            for j in range(nvisit):
                idx = i * nvisit + j

                parinfo[idx]["value"] = float(fit_par["value"][ii])
                parinfo[idx]["fixed"] = _as_bool(fit_par["fixed"][ii])

                if j > 0:
                    parinfo[idx]["tied"] = f"p[{nvisit * i}]"

                if np.isfinite(lo):
                    parinfo[idx]["limited"][0] = True
                    parinfo[idx]["limits"][0] = lo

                if np.isfinite(hi):
                    parinfo[idx]["limited"][1] = True
                    parinfo[idx]["limits"][1] = hi

                params_s.append(float(fit_par["value"][ii]))

            ii += 1

        else:
            # One row per visit.
            for j in range(nvisit):
                idx = i * nvisit + j

                prior = fit_par["prior"][ii]
                p1 = fit_par["p1"][ii]
                p2 = fit_par["p2"][ii]
                lo, hi = prior_bounds(prior, p1, p2, nsigma=5.0)

                parinfo[idx]["value"] = float(fit_par["value"][ii])
                parinfo[idx]["fixed"] = _as_bool(fit_par["fixed"][ii])

                if np.isfinite(lo):
                    parinfo[idx]["limited"][0] = True
                    parinfo[idx]["limits"][0] = lo

                if np.isfinite(hi):
                    parinfo[idx]["limited"][1] = True
                    parinfo[idx]["limits"][1] = hi

                params_s.append(float(fit_par["value"][ii]))
                ii += 1

    return parinfo, np.array(params_s)


def get_initial_walkers(data, theta, meta, fit_par):
    """Generate initial MCMC walker positions from fit_par priors.

    For U priors:
        draw uniformly from [p1, p2]

    For N priors:
        draw from Normal(mean=p1, sigma=p2)

    For X priors:
        raise an error for free parameters, because MCMC needs a finite
        initialization scale.
    """
    nwalkers = meta.run_nwalkers
    ndim = len(theta)

    priors = []
    ii = 0
    nvisit = int(meta.nvisit)

    for i in range(len(fit_par)):
        if ii == len(fit_par):
            break

        fixed = _as_bool(fit_par["fixed"][ii])
        tied = str(fit_par["tied"][ii])

        if not fixed:
            if tied == "-1":
                priors.append(
                    (
                        str(fit_par["prior"][ii]).upper(),
                        float(fit_par["p1"][ii]),
                        float(fit_par["p2"][ii]),
                    )
                )
                ii += 1
            else:
                for _ in range(nvisit):
                    priors.append(
                        (
                            str(fit_par["prior"][ii]).upper(),
                            float(fit_par["p1"][ii]),
                            float(fit_par["p2"][ii]),
                        )
                    )
                    ii += 1
        else:
            ii += 1

    if len(priors) != ndim:
        raise ValueError(
            f"Expected {ndim} priors for MCMC initialization, "
            f"but found {len(priors)}."
        )

    pos = np.zeros((nwalkers, ndim), dtype=float)

    rng = np.random.default_rng()

    for j, (prior, p1, p2) in enumerate(priors):
        if prior == "U":
            if p2 <= p1:
                raise ValueError(
                    f"Uniform prior requires p2 > p1, got p1={p1}, p2={p2}."
                )
            pos[:, j] = rng.uniform(p1, p2, size=nwalkers)

        elif prior == "N":
            if p2 <= 0:
                raise ValueError(
                    f"Normal prior requires sigma p2 > 0, got p2={p2}."
                )
            pos[:, j] = rng.normal(p1, p2, size=nwalkers)

        elif prior == "X":
            raise ValueError(
                "Cannot initialize MCMC walkers for a free parameter with "
                "prior X. Use U or N."
            )

        else:
            raise ValueError(f"Unknown prior type '{prior}'.")

    return pos