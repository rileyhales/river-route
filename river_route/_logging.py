"""One way to build the logger a class writes its progress to, from the log options on a Configs."""

import itertools
import logging
import sys

__all__ = ['PROGRESS', 'build_logger']

PROGRESS = 25
logging.addLevelName(PROGRESS, 'PROGRESS')

# id() is reused after garbage collection, so it cannot name a logger uniquely
_INSTANCE_COUNT = itertools.count()


def build_logger(cfg, kind: str) -> logging.Logger:
    """
    Build a logger for one instance of ``kind`` from the log options on ``configs``.

    Every instance gets its own logger and owns its single handler. The logger does not propagate, because the
    handler is already attached here and propagating to the root logger would print every message twice.

    Args:
        cfg: a Configs, read for log, log_level, log_stream, and log_format
        kind: what is being built, used to name the logger (e.g. 'router', 'network')
    """
    logger = logging.getLogger(f'river_route.{kind}{next(_INSTANCE_COUNT)}')
    logger.propagate = False
    logger.disabled = not cfg.log
    logger.setLevel(cfg.log_level)
    if cfg.log_stream == 'stdout':
        logger.addHandler(logging.StreamHandler(sys.stdout))
    else:
        logger.addHandler(logging.FileHandler(cfg.log_stream))
    logger.handlers[0].setFormatter(logging.Formatter(cfg.log_format))
    return logger
