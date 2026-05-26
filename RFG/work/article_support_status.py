"""
Article support ledger for the first Georgian RFG article.

This file does not add new physics.  It collects the article-facing theorem
exports from the active work files so the draft can be checked from one place.
"""

from p01_core import article_core_theorem
from p02_cosmo import article_cosmology_theorem
from p03_solar import article_solar_theorem
from p04_gw import article_gw_theorem
from p08_cmb import article_cmb_theorem


def first_article_support_status():
    """Return the computational support map for the first article."""
    return {
        "p01_core": article_core_theorem(),
        "p02_cosmo": article_cosmology_theorem(),
        "p03_solar": article_solar_theorem(),
        "p04_gw": article_gw_theorem(),
        "p08_cmb": article_cmb_theorem(),
    }


def first_article_status_summary():
    """Compact status summary for quick terminal checks."""
    status = first_article_support_status()
    return {
        "core": status["p01_core"]["article_status"],
        "cosmology": status["p02_cosmo"]["article_status"],
        "solar": status["p03_solar"]["article_status"],
        "gw": status["p04_gw"]["article_status"],
        "cmb": status["p08_cmb"]["article_status"],
    }


if __name__ == "__main__":
    for section, value in first_article_status_summary().items():
        print(f"\n[{section}]")
        for key, item in value.items():
            print(f"{key}: {item}")
