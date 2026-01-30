# Draft: RNA-seq Platform Enhancements

## Requirements (confirmed)
- Multi-file upload: merge multiple Excel/CSV files into single count matrix
- Tab explanations: educational help for every 15 tabs (biology-friendly)
- Enhanced InterpretationEngine: deeper, context-aware insights

## Technical Decisions
- All Plotly, no matplotlib
- CSS via st.markdown only
- Don't modify backend modules except interpretation_engine.py
- Main app is single file: rnaseq_analysis_platform.py
- Must not break 142 passing tests

## Scope Boundaries
- INCLUDE: Upload UI, merge logic, HELP_TEXTS for all 15 tabs, interpretation_engine.py enhancements
- EXCLUDE: Backend module changes (de_analysis.py, pathway_enrichment.py, etc.), new analysis types

## Open Questions
- Awaiting agent research results for precise code locations
