"""KPF-CC queue subpackage.

Layout:

- :mod:`astroq.queue.kpfcc.queue` — :class:`KPFCC` telescope/instrument
  descriptor consumed by :class:`astroq.splan.SemesterPlanner`,
  :class:`astroq.access.Access`, :class:`astroq.nplan.NightPlanner`, and the
  TTP MILP.
- :mod:`astroq.queue.kpfcc.prep` — data ingestion: pull KPF-CC observing
  blocks (OBs), histories, and allocation info from Keck APIs; format and
  validate the request/custom/allocation frames. Consumed by
  :func:`astroq.driver.kpfcc_prep`.
- :mod:`astroq.queue.kpfcc.starlist` — nightly starlist writer
  (``write_starlist`` / ``format_kpf_row`` / ``pm_correcter``). Consumed
  indirectly via :meth:`KPFCC.write_starlist` from
  :func:`astroq.driver.plan_night`.
"""
