"""Web application subpackage for AstroQ.

The Flask app and all request handlers live in :mod:`astroq.webapp.app`.
Import that module directly (``import astroq.webapp.app as wa``) so the
module-level globals it mutates (``app``, ``uptree_path``, ...) all share
a single namespace.
"""
