"""
Queue subpackage: telescope + instrument descriptors and the factory that
maps ``[global] queue`` -> concrete :class:`~astroq.queue.base.Queue` subclass.

A Queue instance owns all telescope and instrument knowledge for a single
production run (site, slew model, pointing limits, overheads, starlist
writer). It is constructed once by ``SemesterPlanner.__init__`` and reused
by ``Access``, ``NightPlanner``, and the TTP MILP.

Example:
    >>> from configparser import ConfigParser
    >>> import astroq.queue
    >>> cfg = ConfigParser(); cfg.read('config.ini')
    >>> q = astroq.queue.from_config(cfg)
    >>> isinstance(q, astroq.queue.Queue)
    True
"""

from astroq.queue.base import Queue
from astroq.queue.hirescps import HIRESCPS
from astroq.queue.kpfcc import KPFCC

QUEUE_REGISTRY = {
    "hirescps": HIRESCPS,
    "kpfcc": KPFCC,
}


def from_name(name):
    """Build a Queue from its registry name (e.g. ``"hirescps"``).

    Used by the h5 rehydration path where we only have the queue name as a
    string, not a full config object.
    """
    try:
        cls = QUEUE_REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"Unknown queue {name!r}. Known queues: {sorted(QUEUE_REGISTRY)}."
        )
    return cls()


def from_config(config):
    """Build a Queue from a ``[global] queue`` config field."""
    if not config.has_option("global", "queue"):
        raise KeyError(
            "Config is missing [global] queue. Set it to one of "
            f"{sorted(QUEUE_REGISTRY)}."
        )
    return from_name(config.get("global", "queue"))


__all__ = ["Queue", "HIRESCPS", "KPFCC", "QUEUE_REGISTRY", "from_config", "from_name"]
