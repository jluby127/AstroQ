"""Compatibility shim for upstream TTP.

The real telescope/instrument definitions live in :mod:`astroq.queue`. They
are re-exported here so that the vendored TTP source (``astroq.ttp.formatting``,
``astroq.ttp.star``, ``astroq.ttp.model``) and the standalone ``bin/runTTP.py``
keep working with bare imports like ``from astroq.ttp.telescope import Keck1``.

Production code paths in AstroQ build the Queue via
:func:`astroq.queue.from_config` and pass it directly into TTP as the
``observatory`` argument; this shim only matters for one-off scripts and any
remaining bare-import call sites inside the vendored TTP code.
"""

from astroq.queue.hirescps import HIRESCPS as Keck1

# Friendly registry mirroring upstream TTP's ``tel_map``. Add additional
# Queues here as we add more concrete subclasses.
tel_map = {
    "Keck1": Keck1,
}


def create_tel(name):
    """Construct a Queue by upstream-TTP name (``"Keck1"`` ...).

    Returns ``None`` (and prints a warning) for unknown names, matching the
    upstream TTP behavior so existing standalone scripts continue to run.
    """
    cls = tel_map.get(name)
    if cls is None:
        print(
            f"Telescope {name!r} is not in our pre-built classes. See "
            "astroq.queue for the canonical Queue registry."
        )
        return None
    return cls()
