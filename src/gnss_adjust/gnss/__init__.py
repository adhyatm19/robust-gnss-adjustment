"""RINEX-to-baseline integration: RTKLIB invocation, parsing, frames, network assembly.

GNSS carrier-phase processing (ambiguity resolution, orbits, atmosphere) is
delegated to RTKLIB's ``rnx2rtkp``; this subpackage only prepares inputs,
records provenance, parses and validates the resulting baseline vectors and
covariances, and converts them to the repository's network format.
"""
