"""Shared, explicit coordinates: all intervals are zero-based, half-open."""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Transcript:
    id: str
    gene: str
    contig: str
    strand: str
    cds_intervals: tuple[tuple[int, int], ...]
    reference_cds: str


@dataclass(frozen=True)
class Variant:
    id: str
    pos: int
    ref: str
    alleles: tuple[str, str]
    gt: tuple[int, int]
    phased: bool
    phase_set: str | None
    cds_start: int
    cds_end: int
    cds_ref: str
    cds_alleles: tuple[str, str]
    gq: float | None = None
    dp: int | None = None


@dataclass(frozen=True)
class Evidence:
    start: int
    end: int
    status: str
    detail: str = ""


@dataclass
class Locus:
    transcript: Transcript
    sample: str
    variants: list[Variant]
    evidence: list[Evidence]
    warnings: list[str] = field(default_factory=list)
    source: str = ""


class PhaseSpaceError(ValueError):
    """An input cannot be interpreted within the supported model."""
