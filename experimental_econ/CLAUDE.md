# CLAUDE.md — Behavioral Economics PSET (MRes TSE 2025-2026)

## Project Overview

This project synthesizes academic papers and personal notes to produce a structured LaTeX document with original ideas for answering a graduate-level MRes economics problem set on behavioral development economics.

---

## Directory Structure

```
/Users/jorditorresvallverdu/Documents/GitHub/TSE-2025-2026/claude_trial/
├── CLAUDE.md                          ← this file
├── MRES_behavioralDevelopment_PSET.pdf   ← problem set to answer
├── papers/                            ← reference PDFs
│   └── *.pdf
├── my_ideas.rtf                      ← personal notes and ideas
└── output/
    └── pset_ideas.tex                 ← output document (to be generated)
```

---

## Workflow

Execute steps in order. Do not skip ahead. Also, make sure not to parallelize the reading of PDFs or other processes too much. Go one by one. I have an API limit and it will break otherwise. I don't mind if the process is slow, but good. 

### Step 1 — Read the Problem Set

- Open `MRES_behavioralDevelopment_PSET.pdf`
- Identify and list **every question** clearly
- Note the structure: which questions are related, what each one is testing
- Keep all questions in scope — none should be dropped

### Step 2 — Read all papers in `papers/`

For each PDF, extract and internally organize the following:

1. **Behavioral frictions and puzzles** — what anomaly or deviation from standard theory does the paper document?
2. **Behavioral evidence** — in case the paper is empirical, what is the empirical strategy, data, and key result?
3. **Theoretical models** — what model is used or proposed? What are the key mechanisms?
4. **Identification**- what is the main source of identification in the empirical ones - for those that are empirical. 
4. **Research gaps and future lines** — what does the paper leave open, or what does it suggest as next steps?

Cross-reference across papers: identify common themes, complementary findings, and tensions between approaches.

### Step 3 — Read `my_ideas.rtf`

- Understand how Jordi views these papers and what connections he has already drawn
- Use this as the intellectual anchor — the output should reflect and extend these ideas, not replace them
- Note which ideas are most developed vs. rough sketches
- Prioritize sketch of identification and econometric arguments that link RCT variation (and design) to a structural model that derives channels.

### Step 4 — Cross-reference with frontier research

- **Use `WebSearch` and `WebFetch` tools directly** to find recent papers. Do NOT rely on training data or memory for citations — actually search the web. If a subagent is used, verify it performed real web searches.
- Search for recent (2023–2026) NBER working papers and published articles on: information treatments in education, parental beliefs and investments, salience in education decisions, SMS/WhatsApp interventions in schools, behavioral-structural models of education.
- For each paper found, verify it exists by fetching the URL or abstract. **Do not hallucinate citations.**
- Identify where the papers + personal ideas connect to **current NBER working papers** on behavioral and structural economics of education.
- In particular, I am interested in understanding **Information treatments/problems** in education.
- Flag lines of research that appear **genuinely new or underexplored** as of 2025-2026.
- Prioritize novelty: ideas that combine behavioral models with structural approaches are of particular interest.
- **API budget awareness**: This step involves many web calls. Do searches sequentially (not in parallel) and limit to 5–8 targeted searches. Fetch only abstracts/first pages, not full papers. Stop searching once you have 5–10 solid frontier references.

### Step 5 — Generate output

Write `output/pset_ideas.tex` as described below.

---

## Output Specification

**File:** `output/pset_ideas.tex`

**Purpose:** A structured reference document for Jordi to use when writing his own PSET answers. This is NOT a draft of the answers themselves — it is an organized synthesis of key ideas he can develop.

**Structure:**

```latex
\documentclass{article}
% standard packages: amsmath, amssymb, hyperref, geometry, biblatex or natbib
```

One section per PSET question, each containing:

- **Key concepts to mobilize** — relevant behavioral frictions, models, and evidence from the papers
- **Original angles** — new ideas or combinations that go beyond what any single paper does
- **Potential research contributions** — what would be genuinely novel if developed formally
- **Relevant citations** — author (year) inline, no need for full bibliography unless natural

**Tone:** Concise, technical, graduate-level. No padding. Bullet points and short paragraphs preferred over prose.

---

## Constraints

- Do **not** write the PSET answers — only the intellectual raw material to build them
- Do **not** summarize papers exhaustively — focus only on what is relevant to the PSET questions
- Prioritize **novel combinations** of ideas over standard textbook points
- LaTeX must compile cleanly — use only standard packages
- Keep `CLAUDE.md` under 200 lines

---

## Definitions (project-specific)

- **Behavioral frictions**: deviations from standard rationality (present bias, inattention, reference dependence, etc.)
- **Structural approach**: models with explicit primitives, estimated structurally, policy-relevant counterfactuals
- **New line**: an idea not yet mainstream in published journals as of 2025, ideally connectable to an NBER working paper
