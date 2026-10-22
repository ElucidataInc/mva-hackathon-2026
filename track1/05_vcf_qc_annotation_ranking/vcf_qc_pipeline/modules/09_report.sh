#!/usr/bin/env bash
# Module 09 – parse every tool output and render the HTML report.

mod_report() {
    local d="${OUTDIR}/09_report"; mkdir -p "$d"
    "$PYTHON" "${SCRIPTS}/parse_tool_outputs.py" --outdir "$OUTDIR" --out "$d/tool_outputs.json" || return 1
    "$PYTHON" "${SCRIPTS}/build_report.py" --outdir "$OUTDIR" --report "${OUTDIR}/${REPORT_NAME}" --title "$TITLE"
}
