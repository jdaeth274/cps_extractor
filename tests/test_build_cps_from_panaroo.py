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


def test_parse_presence_cell():
    assert mod.parse_presence_cell("") == []
    assert mod.parse_presence_cell("id1;id2") == ["id1", "id2"]
    assert mod.parse_presence_cell('"id1"; "id2"') == ["id1", "id2"]


def test_find_column_flexible_names():
    fields = ["Gene", "Genome Name", "DNA sequence"]
    assert mod.find_column(fields, ["gene_id", "gene"]) == "Gene"
    assert mod.find_column(fields, ["isolate", "genome"]) == "Genome Name"
    assert mod.find_column(fields, ["dna_sequence"]) == "DNA sequence"


def test_read_gene_data_parses_minimal_file(tmp_path):
    gd = tmp_path / "gene_data.csv"
    gd.write_text(
        "Gene,isolate,dna_sequence\n"
        "group_1,isoA,ATGC\n"
        "group_1,isoB,ATGA\n"
        "group_2,isoA,TTTT\n"
    )

    isolates, mapping = mod.read_gene_data(gd)
    assert isolates == {"isoA", "isoB"}
    assert mapping["group_1"]["isoA"] == ["ATGC"]
    assert mapping["group_1"]["isoB"] == ["ATGA"]
