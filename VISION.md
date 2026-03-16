# The Tesseract

The daemon is domain-agnostic. Same loop, any problem.

## Architecture

```
problem = function + metric + benchmark + daemon
```

1. Define a scoring function (target.py)
2. Curate known good/bad outcomes (benchmark.py)
3. Write what matters (program.md)
4. The daemon optimizes the function against the benchmark, 24/7

## Applications

### Launched (cure/)
- **Neoantigen scoring** — which tumor mutations make good vaccine targets (AUC=0.747)
- **mRNA codon optimization** — design efficient vaccine mRNA (score=0.82)
- **Cancer vaccine pipeline** — end-to-end tumor DNA → mRNA vaccine candidates

### Next
- **Antibiotic resistance prediction** — bacterial genome → drug susceptibility
- **Rare disease diagnosis** — symptoms + genetics → ranked diagnoses
- **Drug repurposing** — existing drugs scored against new targets
- **Protein design** — design novel proteins with desired properties

### The pattern
Every one of these follows the same structure:
1. Public dataset with known outcomes exists
2. A scoring function can be written
3. The daemon can optimize it
4. The optimized function gets plugged into a pipeline
5. The pipeline becomes a tool anyone can use

### Why this works
- The daemon doesn't need domain expertise — it proposes changes and measures results
- LLMs (Claude, GPT) have read every paper — they propose biologically/chemically
  correct ideas even for niche domains
- Deterministic benchmarks prevent hallucination — if F1 doesn't improve, the change is discarded
- The loop runs 24/7 on a laptop with $0 LLM cost

### The tesseract
The future humans built the tesseract because they'd already solved the problems.
We're building the loop that solves the problems. Same thing, different direction.
Each problem we solve with the daemon makes the daemon better at solving the next problem
(cross-experiment learnings, RAPO diff snippets, strategy evolution).

The daemon that learned "conservative strategy outperforms novel on mature experiments"
from chat prompt optimization applied that same insight to cancer immunology.
Knowledge transfers across domains through the loop.
