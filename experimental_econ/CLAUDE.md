# CLAUDE.md — Behavioral Economics presentation (MRes TSE 2025-2026)

## Project Overview

This project contains the sources needed to create a presentation in beamer of a 15 min research proposal for a PhD economic course in Experimental Development Economics
---

## Directory Structure

```
/Users/jorditorresvallverdu/Documents/GitHub/TSE-2025-2026/experimental_econ/
├── CLAUDE.md                          ← this file
├── behavioral_pset_torresval.pdf              <- this is the backbone of my research proposal 
└── output/
    └── presentation_rprop_exp.tex                 ← output document (to be generated) + folder
```

---

## Workflow

Execute steps in order. Go one by one. I don't mind if the process is slow, but good.

### Step 1 — behavioral_pset_torresval.pdf

- Open and read `behavioral_pset_torresval.pdf`
- This contains my ideas for the research proposal.
- Make sure to understand the underlying problem, the model proposed and how the data needs to be collected.

### Step 2 — Create a folder "output" inside this directory

### Step 3 - Inside the folder "output" create a presentation for the research proposal

The guidelines to consider are the following, based on the behavioral_pset_torresval.pdf info:

1. **Presentation length** : maximum 15 minutes, so make it simple
2. **Create a beamer**: create a beamer tex file, and consider the following slides and guidelines:
3. **Structure of slides**:
    2.1 Introduction: 1 slide of context, informational problems in education, why is this important (cite relevant papers from structural econometrics: Arcidiacono, Maurel, Aucejo Ransom; Stinebrickner, Wiswall Zafar, Larroucau Rios;) why is this an issue (students seem to be missinformed)
    2.2 Experimental evidence on a concrete issue of information: providing info treatments to parents (bergman, dizon-ross)
    2.3. Caveat: hard to understand the mechanisms behind this last approach: mention the possible mechanisms.
    2.4. Proposal->questions to answer (basically be able to separate mechanisms)
    2.5. The model
    2.6. The RCT
    2.7. Identification based on RCT
    2.8. Power Calculations? (he also asks us to do some of this)
    2.9. Conclusion.

4. **Important guidelines** Be aware that presentation for research must not be crowded in text.
    3.1 This means short sentences better than complete sentences, the slides need to guide the proposal, not the oppsite.
    3.2. Use real citations of papers, don't invent shit. So make sure citations are correct. Main papers are already included in the document behavioral_pset_torresval.
    3.3. Be precise in the math of the model, that is important to explain well.
    3.4. In the RCT, make it visual, maybe with a table, also for identification.

### Step 4 - Output

1. **Output** : save the presentation as "presentation_rprop_exp" inside the output folder, render it into pdf too.

### Step 5- Revise

1. Go through presentation_rprop_exp and revise that guidelines are completed and that format is correct.
2. Propose improvements to clarity and format
3. Save a new document called presentation_rprop_exp_v2 in tex and pdf

END

###### EOF
