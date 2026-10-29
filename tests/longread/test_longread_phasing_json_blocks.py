"""Regression tests for phase blocks that went through JSON and for read alleles that are not calls.

``variant_haplotypes`` maps integer variant positions to alleles. JSON object keys are always text, so a
block written to JSON and read back had ``"1000"`` where ``haplotag_reads`` looked up ``1000``: no read
position matched, every read kept HP=0, and nothing said why. The keys are normalised to integers on input
now, so reloaded blocks tag exactly the reads the originals tag.

An allele that is neither 0 nor 1 (a missing ``allele`` key, a no-call) used to be a vote for haplotype 2 in
``haplotag_reads`` and a reference observation in ``phase_reads``.

All tests run the real functions on small synthetic datasets and on real JSON text; there are no mocks.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import pytest

from metainformant.longread.phasing import haplotyping

POSITIONS = [100, 200, 300, 400, 500, 600]
ANTI_PHASE = [[0, 1, 0, 1, 1, 0], [1, 0, 1, 0, 0, 1]]  # two haplotypes that differ at every variant


def _variants(chrom: str, positions: list[int]) -> list[dict[str, Any]]:
    return [{"chrom": chrom, "position": p, "ref": "A", "alt": "G"} for p in positions]


def _reads(prefix: str, chrom: str, positions: list[int], patterns: list[list[int]], copies: int = 6) -> list[dict]:
    reads = []
    for pattern_index, pattern in enumerate(patterns):
        for copy in range(copies):
            reads.append(
                {
                    "read_id": f"{prefix}{pattern_index}_{copy}",
                    "chrom": chrom,
                    "start": min(positions) - 10,
                    "end": max(positions) + 10,
                    "alleles": [{"position": p, "allele": a} for p, a in zip(positions, pattern)],
                }
            )
    return reads


def _through_json(value: Any) -> Any:
    return json.loads(json.dumps(value))


def _hp_ps(tagged: list[dict[str, Any]]) -> list[tuple[int, int]]:
    return [(r["HP"], r["PS"]) for r in tagged]


class TestBlocksReloadedFromJson:
    def test_json_really_turns_the_position_keys_into_text(self) -> None:
        # the premise of the bug: nothing in the reloaded map can equal an int position
        blocks = haplotyping.phase_reads(_variants("chr1", POSITIONS), _reads("r", "chr1", POSITIONS, ANTI_PHASE))[
            "phase_blocks"
        ]
        reloaded = _through_json(blocks)
        assert all(isinstance(key, int) for key in blocks[0]["variant_haplotypes"])
        assert all(isinstance(key, str) for key in reloaded[0]["variant_haplotypes"])
        assert sorted(map(int, reloaded[0]["variant_haplotypes"])) == POSITIONS

    def test_reloaded_phase_reads_blocks_tag_the_same_reads_as_the_originals(self) -> None:
        reads = _reads("r", "chr1", POSITIONS, ANTI_PHASE)
        blocks = haplotyping.phase_reads(_variants("chr1", POSITIONS), reads)["phase_blocks"]

        original = haplotyping.haplotag_reads(reads, blocks)
        reloaded = haplotyping.haplotag_reads(reads, _through_json(blocks))

        assert {hp for hp, _ in _hp_ps(original)} == {1, 2}  # the originals do tag
        assert _hp_ps(reloaded) == _hp_ps(original)  # and so do the reloaded blocks: they used to tag nothing

    def test_reloaded_build_phase_blocks_blocks_tag_the_same_reads_as_the_originals(self) -> None:
        phased = [
            {"chrom": "chr1", "position": p, "haplotype": h, "quality": 0.9}
            for p, h in zip(POSITIONS, [1, 2, 1, 2, 2, 1])
        ]
        blocks = haplotyping.build_phase_blocks(phased)
        assert len(blocks) == 1
        reads = _reads("r", "chr1", POSITIONS, ANTI_PHASE)

        original = haplotyping.haplotag_reads(reads, blocks)
        reloaded = haplotyping.haplotag_reads(reads, _through_json(blocks))

        assert {hp for hp, _ in _hp_ps(original)} == {1, 2}
        assert _hp_ps(reloaded) == _hp_ps(original)

    def test_blocks_written_to_a_file_and_read_back_keep_their_block_index(self, tmp_path: Path) -> None:
        # two interleaved blocks: PS is the index of the block in the list, which must survive the round trip
        x_positions, y_positions = [100, 300, 500], [200, 400, 600]
        reads = _reads("x", "chr1", x_positions, [[0, 0, 0], [1, 1, 1]]) + _reads(
            "y", "chr1", y_positions, [[0, 1, 0], [1, 0, 1]]
        )
        variants = _variants("chr1", sorted(x_positions + y_positions))
        blocks = haplotyping.phase_reads(variants, reads)["phase_blocks"]
        assert len(blocks) == 2

        path = tmp_path / "blocks.json"
        path.write_text(json.dumps(blocks))
        reloaded = json.loads(path.read_text())

        original = haplotyping.haplotag_reads(reads, blocks)
        again = haplotyping.haplotag_reads(reads, reloaded)
        assert _hp_ps(again) == _hp_ps(original)
        assert {ps for _, ps in _hp_ps(again)} == {0, 1}
        assert all(hp in (1, 2) for hp, _ in _hp_ps(again))

    def test_blocks_on_two_chromosomes_after_json(self) -> None:
        reads = _reads("a", "chr1", POSITIONS, ANTI_PHASE) + _reads("b", "chr2", POSITIONS, ANTI_PHASE)
        variants = _variants("chr1", POSITIONS) + _variants("chr2", POSITIONS)
        blocks = haplotyping.phase_reads(variants, reads)["phase_blocks"]
        assert _hp_ps(haplotyping.haplotag_reads(reads, _through_json(blocks))) == _hp_ps(
            haplotyping.haplotag_reads(reads, blocks)
        )

    def test_the_reloaded_map_is_not_reported_as_missing(self, caplog: pytest.LogCaptureFixture) -> None:
        blocks = haplotyping.phase_reads(_variants("chr1", POSITIONS), _reads("r", "chr1", POSITIONS, ANTI_PHASE))[
            "phase_blocks"
        ]
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            haplotyping.haplotag_reads(_reads("r", "chr1", POSITIONS, ANTI_PHASE), _through_json(blocks))
        assert not [r for r in caplog.records if "variant_haplotypes" in r.getMessage()]

    def test_the_input_blocks_are_not_modified(self) -> None:
        reloaded = _through_json(
            [{"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {"100": 0, "300": 1}}]
        )
        haplotyping.haplotag_reads(_reads("r", "chr1", [100, 300], [[0, 1]]), reloaded)
        assert reloaded[0]["variant_haplotypes"] == {"100": 0, "300": 1}


class TestKeyNormalisation:
    READ = {
        "read_id": "r",
        "chrom": "chr1",
        "start": 50,
        "end": 400,
        "alleles": [{"position": 100, "allele": 0}, {"position": 300, "allele": 1}],
    }

    @staticmethod
    def _block(variant_haplotypes: Any) -> dict[str, Any]:
        return {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": variant_haplotypes}

    def test_text_keys_vote_like_int_keys(self) -> None:
        (as_int,) = haplotyping.haplotag_reads([self.READ], [self._block({100: 0, 300: 1})])
        (as_text,) = haplotyping.haplotag_reads([self.READ], [self._block({"100": 0, "300": 1})])
        assert (as_text["HP"], as_text["PS"]) == (as_int["HP"], as_int["PS"]) == (1, 0)

    def test_text_keys_with_the_opposite_alleles_tag_haplotype_two(self) -> None:
        (tagged,) = haplotyping.haplotag_reads([self.READ], [self._block({"100": 1, "300": 0})])
        assert tagged["HP"] == 2

    def test_int_and_text_keys_can_share_a_block(self) -> None:
        (tagged,) = haplotyping.haplotag_reads([self.READ], [self._block({100: 0, "300": 1})])
        assert tagged["HP"] == 1

    def test_padded_text_and_numpy_integer_keys(self) -> None:
        np = pytest.importorskip("numpy")
        (padded,) = haplotyping.haplotag_reads([self.READ], [self._block({" 100 ": 0, "+300": 1})])
        (numpy_keys,) = haplotyping.haplotag_reads([self.READ], [self._block({np.int64(100): 0, np.int32(300): 1})])
        assert padded["HP"] == numpy_keys["HP"] == 1

    def test_one_position_under_two_spellings_with_one_allele_is_fine(self) -> None:
        (tagged,) = haplotyping.haplotag_reads([self.READ], [self._block({100: 0, "100": 0, "300": 1})])
        assert tagged["HP"] == 1

    def test_one_position_under_two_spellings_with_two_alleles_raises(self) -> None:
        with pytest.raises(ValueError, match="position 100 twice"):
            haplotyping.haplotag_reads([self.READ], [self._block({100: 0, "100": 1})])

    @pytest.mark.parametrize("bad_key", ["chr1:100", "abc", "", "1e2", "100.5", None, 100.5, (100,)])
    def test_a_key_that_is_not_an_integer_position_raises(self, bad_key: Any) -> None:
        with pytest.raises(ValueError, match="not an integer position"):
            haplotyping.haplotag_reads([self.READ], [self._block({bad_key: 0})])

    def test_the_error_names_the_block(self) -> None:
        blocks = [self._block({100: 0}), self._block({"x": 0})]
        with pytest.raises(ValueError, match="phase block 1"):
            haplotyping.haplotag_reads([self.READ], blocks)

    @pytest.mark.parametrize("bad_allele", [2, -1, "1", None, 0.5])
    def test_an_allele_that_is_not_zero_or_one_raises(self, bad_allele: Any) -> None:
        with pytest.raises(ValueError, match="allele is 0 .ref. or 1 .alt."):
            haplotyping.haplotag_reads([self.READ], [self._block({"100": bad_allele})])

    def test_a_map_that_is_not_a_mapping_raises(self) -> None:
        with pytest.raises(TypeError, match="variant_haplotypes must map positions to alleles"):
            haplotyping.haplotag_reads([self.READ], [self._block([100, 300])])

    def test_missing_none_and_empty_maps_still_just_warn(self, caplog: pytest.LogCaptureFixture) -> None:
        blocks = [
            {"chrom": "chr1", "start": 100, "end": 300},
            self._block(None),
            self._block({}),
            self._block({"100": 0, "300": 1}),
        ]
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            (tagged,) = haplotyping.haplotag_reads([self.READ], blocks)
        assert (tagged["HP"], tagged["PS"]) == (1, 3)
        assert any("3 of 4 phase blocks carry no variant_haplotypes" in r.getMessage() for r in caplog.records)


class TestAllelesThatAreNotCalls:
    BLOCK = {"chrom": "chr1", "start": 100, "end": 300, "variant_haplotypes": {100: 0, 300: 0}}

    @staticmethod
    def _read(alleles: list[dict[str, Any]]) -> dict[str, Any]:
        return {"read_id": "r", "chrom": "chr1", "start": 50, "end": 400, "alleles": alleles}

    @pytest.mark.parametrize(
        "uncalled", [{"position": 100}, {"position": 100, "allele": None}, {"position": 100, "allele": 2}]
    )
    def test_haplotag_does_not_count_a_missing_or_no_call_allele_as_a_vote_for_haplotype_two(
        self, uncalled: dict[str, Any]
    ) -> None:
        # one real vote for haplotype 1 and one uncalled observation: it used to be a 1-1 tie, HP 0
        (tagged,) = haplotyping.haplotag_reads([self._read([uncalled, {"position": 300, "allele": 0}])], [self.BLOCK])
        assert (tagged["HP"], tagged["PS"]) == (1, 0)

    def test_haplotag_reports_how_many_alleles_it_did_not_count(self, caplog: pytest.LogCaptureFixture) -> None:
        reads = [self._read([{"position": 100}, {"position": 300, "allele": 7}]) for _ in range(3)]
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            tagged = haplotyping.haplotag_reads(reads, [self.BLOCK])
        assert [r["HP"] for r in tagged] == [0, 0, 0]
        assert any("6 read alleles" in r.getMessage() for r in caplog.records)

    def test_phase_reads_does_not_turn_a_missing_allele_into_a_reference_observation(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        variants = _variants("chr1", [100, 200])
        # allele-less observations used to read as allele 0 and link the two variants as in phase
        reads = [
            {
                "read_id": f"r{i}",
                "chrom": "chr1",
                "start": 50,
                "end": 250,
                "alleles": [{"position": 100}, {"position": 200}],
            }
            for i in range(5)
        ]
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            result = haplotyping.phase_reads(variants, reads)
        assert result["haplotype_assignments"] == {} and result["phase_blocks"] == []
        assert result["n_phased_reads"] == 0
        assert any("10 read alleles" in r.getMessage() for r in caplog.records)

    def test_phase_reads_still_uses_the_calls_next_to_an_uncalled_allele(self) -> None:
        variants = _variants("chr1", [100, 200, 300])
        reads = [
            {
                "read_id": f"r{i}",
                "chrom": "chr1",
                "start": 50,
                "end": 350,
                "alleles": [{"position": 100, "allele": 0}, {"position": 200}, {"position": 300, "allele": 0}],
            }
            for i in range(4)
        ]
        result = haplotyping.phase_reads(variants, reads)
        assert result["n_phased_reads"] == 4
        (block,) = result["phase_blocks"]
        assert block["n_variants"] == 2 and sorted(block["variant_haplotypes"]) == [100, 300]
