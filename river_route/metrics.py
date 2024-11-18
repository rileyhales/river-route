import numpy as np


def mean_error(y_true, y_pred):
    return np.mean(y_true - y_pred)


def mean_absolute_error(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))


def mean_square_error(y_true, y_pred):
    return np.mean((y_true - y_pred) ** 2)


def pearson_r(y_true, y_pred):
    return np.corrcoef(y_true, y_pred)[0, 1]


def kling_gupta_efficiency_2012(y_true, y_pred):
    pearson_r = np.corrcoef(y_true, y_pred)[0, 1]
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    std_true = np.std(y_true)
    std_pred = np.std(y_pred)

    beta = mean_pred / mean_true
    gamma = (mean_pred / std_pred) / (mean_true / std_true)

    term1 = np.power(pearson_r - 1, 2)
    term2 = np.power(beta - 1, 2)
    term3 = np.power(gamma - 1, 2)
    return 1 - np.sqrt(term1 + term2 + term3)


def me(y_true, y_pred):
    return mean_error(y_true, y_pred)


def mae(y_true, y_pred):
    return mean_absolute_error(y_true, y_pred)


def mse(y_true, y_pred):
    return mean_square_error(y_true, y_pred)


def r(y_true, y_pred):
    return pearson_r(y_true, y_pred)


def kge2012(y_true, y_pred):
    return kling_gupta_efficiency_2012(y_true, y_pred)
