"""Operational and developer scripts shipped with AstroQ.

Each module here is runnable via ``python -m astroq.scripts.<name>``:

- ``check_night_plans`` — regression check that compares per-band night-plan
  outputs against the script files. Invoked by the production Makefiles.
- ``ttp_keck1`` — standalone TTP runner for HIRES-CPS on Keck-I. Debug /
  planning-meeting tool.
- ``update_backup_OBs`` — update backup OB JSON files with keyword overrides.

These are not part of the library API. Internal tools only.
"""
