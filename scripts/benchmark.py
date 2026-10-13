"""A small reproducible engine-only benchmark; not a real-genome performance claim."""
from statistics import median
from time import perf_counter

from phasespace.analysis import analyze
from phasespace.models import Locus, Transcript, Variant

sequence = "ATG" + "GCT" * 998 + "TAA"
variants = []
for i in range(8):
    position = 10 + i * 90
    reference = sequence[position]
    alternate = "A" if reference != "A" else "C"
    variants.append(Variant(str(i), position, reference, (reference, alternate),
                            (0, 1), False, None, position, position + 1,
                            reference, (reference, alternate), 60, 30))
locus = Locus(Transcript("BENCH.1", "BENCH", "synthetic", "+", ((0, len(sequence)),), sequence),
              "SYNTHETIC", variants, [])
durations = []
for _ in range(5):
    started = perf_counter()
    result = analyze(locus)
    durations.append(perf_counter() - started)
print(f"{len(sequence)} coding bases, 8 blocks, {result['scenario_count']} arrangements")
print(f"Median engine time: {median(durations):.4f} s; range {min(durations):.4f}–{max(durations):.4f} s")
print("Synthetic engine-only benchmark; excludes input parsing and report generation.")
