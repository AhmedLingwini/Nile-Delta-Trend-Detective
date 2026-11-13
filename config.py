"""Shared settings for the Nile Delta Trend Detective scripts."""

# Seasons as 0-based month indices (0 = January ... 11 = December).
# DJF uses January, February and December of the SAME calendar year.
SEASONS = {
    "annual": list(range(12)),
    "DJF": [11, 0, 1],
    "MAM": [2, 3, 4],
    "JJA": [5, 6, 7],
    "SON": [8, 9, 10],
}
