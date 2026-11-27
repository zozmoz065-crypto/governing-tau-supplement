# Governing τ — replication materials

Code and data for the manuscript **“Architecting Organizational Cyber Resilience: A Governance-Centric Framework for Threat-Intelligence Absorption and Policy Translation Latency”** (Mohammed & Thakur, Shoolini University).

The online appendix is `Online_Appendix.pdf`. No organizational data are used. The code regenerates the numerical results in Section 5.6 and the scoping counts in Table 1B.

## Contents

| File | Reproduces |
|---|---|
| `Online_Appendix.pdf` | Online Appendix A (governance-maturity instrument), B (worked τ example) and C (derivations of Results 1b, 3 and 5) |
| `sim.py` | Discrete-event simulation of the triage queue (Results 1–2); defines the simulator used by `figs2.py` |
| `figs.py` | Full pipeline: 10-replication validation (Table 4A), Figures 3–5 and the Result 1b grid check; writes `figs.json` |
| `figs2.py` | Fast path: Figures 4–5 and the grid check, reusing the stored `figs.json` |
| `fig2.py` | Figure 2 (causal and queueing structure) |
| `calib.py` | Calibration against public KEV exposure evidence (Table 4B) → `calib.json` |
| `figs.json` | Stored simulation outputs used for Table 4A (10 replications × 200,000 arrivals per utilization level) |
| `openalex_queries.txt` | Query URLs for the OpenAlex scoping check (Table 1B), run on 27 November 2026 |

## How to run

```bash
pip install -r requirements.txt
python sim.py      # quick validation run (prints analytic vs simulated values)
python figs.py     # full run: Table 4A, Figures 3–5, grid check (several minutes)
python figs2.py    # fast: Figures 4–5 and grid check from stored figs.json
python fig2.py     # Figure 2
python calib.py    # Table 4B
```

The random seed is fixed (`20261127`). `figs.py` reruns the 10-replication validation; `figs2.py` reuses the stored results.

## Re-running the scoping check

Open each URL in `openalex_queries.txt`; the record count is the `meta.count` field of the JSON response. Counts grow as OpenAlex indexes new records.

## Licence

MIT — see `LICENSE`.
