# 06. Hackathon rules and organiser clarifications

A summary of the official rules, the submission forms and scorer in the Space source, the methods template, and the organisers' answers in the Community tab. Based on Space commit `aeeef5ad`. The official Rules take precedence over this summary.

## Summary

- Both tracks close 24 Oct 2026, 23:59 UTC.
- Track 1: up to 6 CSV uploads per Hugging Face user, each with a GitHub URL and a report. The leaderboard shows each display name's best score.
- Track 2: up to 3 uploads per team from one member; only the latest is judged. Each needs a GitHub URL, a 3-minute pitch video and a report. Approved drugs only; combinations allowed.
- The methods description has one required field: the generative-AI provider, tier and data-handling setting. Use a business or API tier for patient-specific data.
- Repositories must be public after close. Incidental findings stay out of public repositories while thread #20 is open.
- Delete all data by 23 Nov 2026.
- The organisers will not confirm phase (#24). Only FASTQ is provided, no BAM/CRAM (#25). No clinical data beyond HPO terms (#14).

## 1. Timeline

| Item | Date |
|---|---|
| Launch and data release | 24 Aug 2026 |
| Submissions open | 25 Aug 2026 |
| Submissions close, leaderboard frozen | 24 Oct 2026, 23:59 UTC |
| Judging | 24 Oct to mid-December 2026 |
| Winners announced | 18 Dec 2026 (moved from 25 Nov, thread #26) |
| Data deletion deadline | 23 Nov 2026 (30 days after close) |

The data are shared under WCG IRB protocol #20252010. Participants must be 18 or older and register individually. How the two tracks combine into placement is not yet announced (#11).

## 2. Submissions

### 2.1 Track 1: variant prediction

| Field | Required | Check |
|---|---|---|
| Team / display name | No (defaults to username) | Groups leaderboard entries; teammates must type it identically (#19) |
| GitHub URL | Yes | Starts with `https://github.com/` |
| Predictions | Yes | `.csv` |
| Report | Yes | `.pdf` or `.md` |

CSV header (from `track1_submission_template.csv`):

```
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

- `proband_id` must be `PROBAND01`. At most 10 rows. `epcr` in (0, 1]. `finding_type` is `primary` (default) or `secondary`.
- A compound heterozygous pair goes on one row, in the `_1` and `_2` columns. The scorer compares each row's variants with the answer key as a set.
- Coordinates are GRCh38. The template uses `chr15`; the VCF has no prefix (#17). The scorer strips whitespace but does not normalise the prefix, so follow the template.
- REF and ALT are uppercased. Indels must match the key's representation (left alignment, anchor base).
- Name files `<username-or-team>_<model>.csv` and `<username-or-team>_track1_report.pdf` (or `.md`).

Scoring (`evaluation.py`):

| Rank of full match | Points |
|---|---|
| 1 | 100 |
| 2-3 | 50 |
| 4-5 | 25 |
| 6-10 | 10 |
| One allele of the pair only | half |

Rows are sorted by EPCR, ties by file order. F-max is computed per variant across every EPCR threshold, so a secondary row at or above the primary row's EPCR pushes the primary down and adds false positives. The methods write-up is judged separately after close and "counts just as much" (#18).

The repository must hold documented, reproducible code and configuration. Submission files are optional; a methods description is recommended (PDF is fine, #10). Organisers may rerun submissions.

### 2.2 Track 2: drug repurposing

| Field | Required | Check |
|---|---|---|
| Team / display name | No | Used for the announcement |
| GitHub URL | Yes | Starts with `https://github.com/` |
| Pitch video URL | Yes | 3 minutes, YouTube or Vimeo |
| Report | Yes | `.pdf` or `.md` |
| Notes for judges | No | Free text |

- One member submits; only the latest of three attempts is judged (#7, #10).
- The report must characterise the variant mechanism (loss or gain of function, pathway, downstream effect).
- Only existing, market-approved drugs, alone or in explained combinations. Investigational compounds are out (#5). The regulator defining "approved" is not specified.
- Weights: Scientific Rigor 35, Potential Impact 25, Innovation 25, Scalability 15. A separate Innovation/Community award covers "creativity, community impact, or patient-centered design".
- The About tab calls submissions "hypotheses for follow-up, not evidence that a medicine works".

### 2.3 Methods description

`methods_description_form.xlsx` has one sheet per track. Only the generative-AI row is mandatory, but a missing description may be scored lower (Innovation and Scalability in Track 1, Rigor and Innovation in Track 2). The upload forms do not carry it, so it goes in the report or the repository.

Track 1 (one per model, up to 6): team name; model number; approach; **generative AI provider, tier and data setting**; automated or manually curated output (automated preferred); details of any manual review; public or proprietary data (public preferred); public data used; proprietary data used; whether compound heterozygous pairs can be output; handling of secondary findings; run time and cost; abstract of up to 500 words with strengths and limitations.

Track 2 (one per team): team name; approach from variant to drug; **generative AI provider, tier and data setting**; automated or manual candidate identification; details of manual curation; public or proprietary sources (public preferred); public sources used; proprietary sources used; how the mechanism was characterised; time and effort; abstract of up to 500 words.

## 3. Organiser clarifications

From Community threads #1-#25. Threads about access bugs, templates, networking and off-topic posts are left out.

| Thread | Topic | Answer |
|---|---|---|
| #2 | Third-party LLMs | A service that returns results without taking rights over the data is a processor; one that can use, keep or learn from it is a recipient. Allowed if there is no training on inputs or outputs and retention is limited (zero retention not required). Record provider and configuration. 25 Sep: use a business or API tier for patient-specific data "even with training disabled" |
| #4 | Repository | Private during the event, public after close under CC BY 4.0. Must hold reproducible code and configuration |
| #5 | Combinations, investigational drugs | Combinations of approved drugs allowed with reasons for each drug and the combination. Investigational drugs out |
| #7 | Track 2 resubmission | Quota raised from 1 to 3; the latest counts |
| #10 | Announcements | Generative-AI question added to the template; an updated report can be committed without resubmitting; repositories public "once final evaluations begin"; PDF methods description accepted |
| #11 | Prizes across tracks | Still being decided |
| #13 | Dataset citation | Needed only in later publications, not in reports |
| #14 | Extra clinical data | None; consent covers genomic data and HPO terms only |
| #15 | Deletion confirmation | The Rules address accepts mail only from Synapse-registered addresses |
| #16 | Institution | Not required |
| #17 | Reference genome | GRCh38; no MD5 available; reference-agnostic tools not required. A participant read from the VCF header: no-alt + hs38d1, no `chr` prefix, Sentieon 202308.02, INFO/AF as within-sample fraction |
| #18 | Leaderboard | Many perfect scores; the write-up "counts just as much" |
| #19 | Duplicate entries | Caused by mismatched team names; organisers removed the duplicate and restored the quota |
| #20 | Incidental findings | Referred to privacy and IRB leads. Meanwhile keep such repositories private and do not contact the family. Open |
| #21 | Consumer LLM plans | "May be acceptable" if the terms meet #2; otherwise use an API tier. The 25 Sep update in #2 is stricter |
| #22 | External phasing data | Any suitable dataset if described; All of Us suggested as a free alternative to UK Biobank |
| #24 | Phase | Not confirmed before close, even privately. Public information may be used |
| #25 | Reads | FASTQ only; align locally to GRCh38. Illumina DNA PCR-Free on NovaSeq. Insert size to follow |
| #26 | Announcement date | Winners will be announced on 18 Dec 2026. The submission deadline is unchanged |

## 4. Data handling

- Delete all copies within 30 days of close (by 23 Nov 2026): local machines, cloud instances, notebooks, private repositories and derived data.
- This covers raw files and any subset or reformat (including indexes), genotype-scale intermediate tables, caches and notebook state, local prompt logs with pasted variant data, and any model trained on the genome (#2).
- You may keep the submission CSV, HPO terms, gene and pathway rankings, mechanism analyses, drug candidates, code, reports and the pitch.
- Confirm deletion by email. The Rules give `RarediseaserealkidMVAhackathon2026@synapse.org`; the access form gives `MVAHackathon2026@synapse.org`. Send to both from a Synapse-registered address (#15).
- Do not reshare the data. Everyone who handles it must register. Report suspected disclosure through the Sage Help Center (sagebionetworks.jira.com/servicedesk/customer/portal/20). Do not contact the family or MVA Society contacts.

## 5. Publication and licence

- Submissions (predictions, code, reports) are released under CC BY 4.0 with participant attribution.
- No manuscript using the dataset may go to a peer-reviewed venue until the organisers publish their summary. Code and derived outputs can be shared any time. Abstracts and posters need written approval.
- Any public communication must carry the acknowledgement paragraph from the Rules.
- Nothing may re-identify the child or family.
- Cite the dataset in later publications using the reference on the Synapse page (syn76251147).

## Open questions

1. Which deletion address is correct.
2. How the tracks combine into prizes (#11).
3. Incidental findings in public repositories (#20).
4. The `chr` convention in the answer key.
5. Reference MD5 and insert size (#17, #25).
6. Which regulator defines "market-approved". Listing FDA and EMA status per drug avoids the question.

## References

1. Hackathon Space `tabs/rules.py` (Official Rules). https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026 (commit aeeef5ad).
2. Same Space: `tabs/about.py`, `tabs/submit_track1.py`, `tabs/submit_track2.py`, `tabs/faq.py`, `utils.py`, `config.py`, `evaluation.py` (adapted from Stenton et al. 2024, CAGI6), `static/templates/track1_submission_template.csv`, `static/templates/methods_description_form.xlsx`.
3. Space Community discussions #1-#25. https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions.
4. Synapse project syn76251147. https://www.synapse.org/Synapse:syn76251147.
5. Dataset card, `SageBio/mva-hackathon-2026-data`. https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data.
