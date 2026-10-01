import numpy as np
from numpy.typing import ArrayLike

__all__ = [
    'mean_error',
    'mean_absolute_error',
    'mean_square_error',
    'pearson_r',
    'kling_gupta_efficiency_2012',
    'me',
    'mae',
    'mse',
    'kge2012',
]


def mean_error(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """The mean of the simulated minus the observed values: positive when the simulation is too high on average."""
    return float(np.mean(np.asarray(y_pred) - np.asarray(y_true)))


def mean_absolute_error(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """The mean of the absolute differences between the simulated and the observed values."""
    return float(np.mean(np.abs(np.asarray(y_pred) - np.asarray(y_true))))


def mean_square_error(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """The mean of the squared differences between the simulated and the observed values."""
    return float(np.mean((np.asarray(y_pred) - np.asarray(y_true)) ** 2))


def pearson_r(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """The Pearson correlation coefficient of the observed and the simulated values."""
    return float(np.corrcoef(np.asarray(y_true), np.asarray(y_pred))[0, 1])


def kling_gupta_efficiency_2012(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """
    Kling-Gupta efficiency per Kling et al. (2012): 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2), with r the
    correlation, beta the ratio of the simulated to the observed mean, and gamma the ratio of the simulated to the
    observed coefficient of variation. 1 is a perfect simulation. NaN when either series is constant or has a mean of
    zero, where a coefficient of variation is undefined.
    """
    observed, simulated = np.asarray(y_true), np.asarray(y_pred)
    mean_observed, mean_simulated = observed.mean(), simulated.mean()
    std_observed, std_simulated = observed.std(), simulated.std()
    if std_observed == 0 or std_simulated == 0 or mean_observed == 0 or mean_simulated == 0:
        return float('nan')
    r = pearson_r(observed, simulated)
    beta = mean_simulated / mean_observed
    gamma = (std_simulated / mean_simulated) / (std_observed / mean_observed)
    return float(1 - np.sqrt((r - 1) ** 2 + (beta - 1) ** 2 + (gamma - 1) ** 2))


me = mean_error
mae = mean_absolute_error
mse = mean_square_error
kge2012 = kling_gupta_efficiency_2012
