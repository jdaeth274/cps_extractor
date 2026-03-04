import importlib.util
from pathlib import Path

from Bio.SeqFeature import SeqFeature, FeatureLocation


spec = importlib.util.spec_from_file_location(
    "build_cps_from_panaroo", Path("bin/build_cps_from_panaroo.py")
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_is_excluded_pseudo():
    feat = SeqFeature(FeatureLocation(0, 10), type="CDS", qualifiers={"pseudo": [""]})
    assert mod.is_excluded(feat) is True


def test_is_excluded_transposon():
    feat = SeqFeature(
        FeatureLocation(0, 10),
        type="CDS",
        qualifiers={"product": ["IS transposase"]},
    )
    assert mod.is_excluded(feat) is True


def test_is_not_excluded_regular_cps_gene():
    feat = SeqFeature(
        FeatureLocation(0, 10),
        type="CDS",
        qualifiers={"gene": ["wciP"], "product": ["capsular polysaccharide biosynthesis protein"]},
    )
    assert mod.is_excluded(feat) is False


def test_sanitize_isolate_name():
    assert mod.sanitize_isolate_name("iso 1/a") == "iso_1_a"
