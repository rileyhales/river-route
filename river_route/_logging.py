import logging
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .configs.Configs import Configs

__all__ = ['PROGRESS', 'configure_logging']

PROGRESS = 25
logging.addLevelName(PROGRESS, 'PROGRESS')

_PACKAGE_LOGGER = logging.getLogger('river_route')
_HANDLER_NAME = 'river_route'  # the name of the one handler configure_logging adds, so it replaces only its own


def configure_logging(configs: Configs) -> None:
    """
    Point the river_route logger at the log options of ``configs``: one handler writing ``log_format`` to
    ``log_stream`` at ``log_level``, which replaces the handler of the previous call, or no output when ``log`` is off.
    Every module logs to a child of this logger, so the options apply to all of them. They apply to the whole process,
    so the Router built last sets them. The logger does not propagate, so a root logger the application configured
    does not print each line twice.

    Args:
        configs: a Configs, read for log, log_level, log_stream, and log_format
    """
    for handler in [handler for handler in _PACKAGE_LOGGER.handlers if handler.get_name() == _HANDLER_NAME]:
        _PACKAGE_LOGGER.removeHandler(handler)
        handler.close()
    _PACKAGE_LOGGER.propagate = False
    if not configs.log:
        _PACKAGE_LOGGER.setLevel(logging.CRITICAL + 1)
        return
    _PACKAGE_LOGGER.setLevel(configs.log_level)
    stream = configs.log_stream
    handler = logging.StreamHandler(sys.stdout) if stream == 'stdout' else logging.FileHandler(stream)
    handler.set_name(_HANDLER_NAME)
    handler.setFormatter(logging.Formatter(configs.log_format))
    _PACKAGE_LOGGER.addHandler(handler)
    return
