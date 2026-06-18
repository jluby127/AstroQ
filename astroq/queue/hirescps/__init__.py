"""HIRES-CPS queue subpackage.

Layout:

- :mod:`astroq.queue.hirescps.queue` — :class:`HIRESCPS` telescope/instrument
  descriptor consumed by :class:`astroq.splan.SemesterPlanner`,
  :class:`astroq.access.Access`, :class:`astroq.nplan.NightPlanner`, and the
  TTP MILP.
- :mod:`astroq.queue.hirescps.prep` — data ingestion: pull HIRES-CPS request
  sheets and JUMP past-history; build the allocation CSV from the Keck
  schedule. Consumed by :func:`astroq.driver.hirescps_prep`.
- :mod:`astroq.queue.hirescps.starlist` — nightly starlist writer
  (``write_starlist`` / ``format_hires_row``). Consumed indirectly via
  :meth:`HIRESCPS.write_starlist` from :func:`astroq.driver.plan_night`.
"""
