import numpy as np

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


def mean_error(y_true, y_pred):
    return np.mean(y_true - y_pred)


def mean_absolute_error(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))


def mean_square_error(y_true, y_pred):
    return np.mean((y_true - y_pred) ** 2)


def pearson_r(y_true, y_pred):
    return np.corrcoef(y_true, y_pred)[0, 1]


def me(y_true, y_pred):
    return mean_error(y_true, y_pred)


def mae(y_true, y_pred):
    return mean_absolute_error(y_true, y_pred)


def mse(y_true, y_pred):
    return mean_square_error(y_true, y_pred)


def kling_gupta_efficiency_2012(y_true, y_pred):
    """
    Kling-Gupta efficiency per Kling et al. (2012): 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2), with r the
    correlation, beta the ratio of the simulated to the observed mean, and gamma the ratio of the simulated to the
    observed coefficient of variation. 1 is a perfect simulation. NaN when either series is constant or has a mean of
    zero, where a coefficient of variation is undefined.
    """
    pr = pearson_r(y_true, y_pred)
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    std_true = np.std(y_true)
    std_pred = np.std(y_pred)

    if std_true == 0 or std_pred == 0 or mean_true == 0 or mean_pred == 0:
        return np.nan
    beta = mean_pred / mean_true
    gamma = (std_pred / mean_pred) / (std_true / mean_true)

    term1 = np.power(pr - 1, 2)
    term2 = np.power(beta - 1, 2)
    term3 = np.power(gamma - 1, 2)
    return 1 - np.sqrt(term1 + term2 + term3)


def kge2012(y_true, y_pred):
    return kling_gupta_efficiency_2012(y_true, y_pred)
