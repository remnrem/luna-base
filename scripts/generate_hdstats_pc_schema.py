#!/usr/bin/env python3
"""Generate the HDSTATS-PC registry header and reusable V1 JSON recipe.

The editable source is docs/hdstats/hdstats.pc.  Ordinary three-column rows
declare the feature panels; COMPONENT rows preserve the legacy projection
registry and determine each panel's default V1 component count. The V1 recipe
does not assign fitted component identities or orientations. Linked features
are selected from model loadings at projection time.
"""

import argparse
import json
import re
from pathlib import Path


def cpp(value: str) -> str:
    return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'


def feature_key(variable: str, stratum: str) -> str:
    if not stratum or stratum == ".":
        return variable
    parts = []
    for field in stratum.split(";"):
        factor, level = field.split("/", 1)
        parts.append(f"{factor}_{level}")
    return variable + "." + ".".join(parts)


METRIC_BRIEFS = {
    "ACF": "Autocorrelation",
    "BETWEEN_SD": "Between-block standard deviation",
    "BOUT_TRANS_DENS": "Bout-transition density",
    "BOUT_TRANS_DUR": "Bout-transition duration",
    "BOUT_TRANS_MIN_CONF": "Bout-transition minimum confidence",
    "BOUT_TRANS_N": "Bout-transition count",
    "BOUT_TRANS_PEAK_H": "Bout-transition peak entropy",
    "BOUT_TRANS_TOT_PS": "Bout-transition total posterior shift",
    "COM": "Centre of mass",
    "CONF": "Posterior confidence",
    "CONF_DIFF_UN_ST": "Confidence difference, unstable versus stable",
    "CORR_H_PS": "Entropy/posterior-shift correlation",
    "FRAC_LOW_CONF": "Fraction low-confidence",
    "H": "Hypnodensity entropy",
    "H_P90": "90th-percentile hypnodensity entropy",
    "H_RATIO_UN_ST": "Entropy ratio, unstable versus stable",
    "HD_DENS": "Expected state-set ambiguity density",
    "HD_MINS": "Expected jointly ambiguous minutes",
    "HD_PCT": "Hypnodensity percentage",
    "HD_REM_LAT": "REM latency",
    "HD_REM_LAT2": "REM latency II",
    "HD_SME": "Sleep-maintenance efficiency",
    "HD_SOL": "Sleep-onset latency",
    "HD_SPT": "Sleep-period time",
    "HD_WASO": "Wake after sleep onset",
    "MG": "Posterior probability margin",
    "PCT_STABLE": "Percentage stable",
    "PS": "Posterior shift",
    "PS_LAG": "Posterior-shift lag",
    "PS_P90": "90th-percentile posterior shift",
    "PS_RATIO_UN_ST": "Posterior-shift ratio, unstable versus stable",
    "P_N2": "N2 posterior probability",
    "P_R": "REM posterior probability",
    "P_W": "Wake posterior probability",
    "SIM": "Posterior similarity",
    "SLOPE": "Overnight slope",
    "T50": "Median time",
    "TRANS_DENS": "Transition density",
    "TRANS_DUR": "Transition duration",
    "TRANS_MIN_CONF": "Transition minimum confidence",
    "TRANS_N": "Transition count",
    "TRANS_N_MATCHED": "Matched transition count",
    "TRANS_PCT_MATCHED": "Percentage matched transitions",
    "TRANS_PEAK_H": "Transition peak entropy",
    "TRANS_TOT_PS": "Transition total posterior shift",
    "U_COVERAGE": "Ultradian coverage",
    "U_MAG": "Ultradian magnitude",
    "U_MEAN": "Ultradian mean",
    "U_MOD": "Ultradian modulation",
    "U_MOD_COS": "Cosine ultradian modulation",
    "U_MOD_SIN": "Sine ultradian modulation",
    "U_NCYC": "Ultradian cycle count",
    "U_PARTIAL_R2": "Partial ultradian R-squared",
    "U_R2": "Ultradian R-squared",
    "U_SLOPE": "Ultradian slope",
    "WITHIN_SD": "Within-block standard deviation",
}


METRIC_LABELS = {
    "ACF": "ACF",
    "BETWEEN_SD": "between-block SD",
    "BOUT_TRANS_DENS": "bout trans density",
    "BOUT_TRANS_DUR": "bout trans duration",
    "BOUT_TRANS_MIN_CONF": "bout trans min confidence",
    "BOUT_TRANS_N": "bout trans count",
    "BOUT_TRANS_PEAK_H": "bout trans peak entropy",
    "BOUT_TRANS_TOT_PS": "bout trans total PS",
    "COM": "CoM",
    "CONF": "Confidence",
    "CONF_DIFF_UN_ST": "Confidence diff (U/S)",
    "CORR_H_PS": "Entropy-PS corr",
    "FRAC_LOW_CONF": "Low confidence",
    "H": "Entropy",
    "H_P90": "Entropy p90",
    "H_RATIO_UN_ST": "Entropy ratio (U/S)",
    "HD_DENS": "ambiguity density",
    "HD_MINS": "ambiguous minutes",
    "HD_PCT": "HD %",
    "HD_REM_LAT": "REM latency",
    "HD_REM_LAT2": "REM latency II",
    "HD_SME": "Sleep-maintenance eff.",
    "HD_SOL": "Sleep-onset latency",
    "HD_SPT": "Sleep-period time",
    "HD_WASO": "WASO",
    "MG": "Margin",
    "PCT_STABLE": "Stable %",
    "PS": "PS",
    "PS_LAG": "PS lag",
    "PS_P90": "PS p90",
    "PS_RATIO_UN_ST": "PS ratio (U/S)",
    "P_N2": "N2 probability",
    "P_R": "REM probability",
    "P_W": "Wake probability",
    "SIM": "Similarity",
    "SLOPE": "Slope",
    "T50": "T50",
    "TRANS_DENS": "transition density",
    "TRANS_DUR": "transition duration",
    "TRANS_MIN_CONF": "transition min confidence",
    "TRANS_N": "transition count",
    "TRANS_N_MATCHED": "matched transitions",
    "TRANS_PCT_MATCHED": "matched transitions (%)",
    "TRANS_PEAK_H": "transition peak entropy",
    "TRANS_TOT_PS": "transition total PS",
    "U_COVERAGE": "ultradian coverage",
    "U_MAG": "Ultradian magnitude",
    "U_MEAN": "ultradian mean",
    "U_MOD": "Ultradian mod",
    "U_MOD_COS": "ultradian mod cos",
    "U_MOD_SIN": "ultradian mod sin",
    "U_NCYC": "ultradian cycle count",
    "U_PARTIAL_R2": "partial ultradian R2",
    "U_R2": "Ultradian R2",
    "U_SLOPE": "ultradian slope",
    "WITHIN_SD": "within-block SD",
}


def duration_label(seconds: str, long: bool = False) -> str:
    value = int(seconds)
    if value % 3600 == 0:
        number = value // 3600
        return f"{number} hour" + ("s" if long and number != 1 else "") if long else f"{number * 60}min"
    number = value // 60
    return f"{number} minutes" if long else f"{number}min"


def stratum_brief(stratum: str) -> str:
    if not stratum or stratum == ".":
        return ""
    labels: list[str] = []
    for field in stratum.split(";"):
        factor, level = field.split("/", 1)
        if factor == "NSS":
            labels.append(f"{level}-stage")
        elif factor == "SS":
            labels.append(f"scored {level}")
        elif factor == "MIX":
            labels.append(f"{level.replace('_', '/') } mixture")
        elif factor in {"TRANS", "BOUT_TRANS"}:
            labels.append(f"{factor.lower().replace('_', ' ')} {level.replace('to', ' to ')}")
        elif factor == "STABLE":
            labels.append("stable" if level == "1" else "unstable")
        elif factor == "LRLAG":
            labels.append(f"long-range lag {duration_label(level, long=True)}")
        elif factor == "LRWIN":
            labels.append(f"long-range window {duration_label(level, long=True)}")
        elif factor == "VAR":
            labels.append(f"{level.replace('_', ' ')} variability")
        elif factor == "UQ":
            labels.append(f"ultradian {level}")
        elif factor == "UMETRIC":
            labels.append(f"ultradian {level}")
        elif factor == "UVAR":
            labels.append(f"ultradian {level} variability")
        else:
            labels.append(f"{factor.lower()} {level}")
    return "; ".join(labels)


def feature_brief(variable: str, stratum: str) -> str:
    factors = dict(field.split("/", 1) for field in stratum.split(";") if stratum and stratum != ".")
    if variable in {"HD_DENS", "HD_MINS"} and "SS" in factors:
        stage = factors["SS"]
        return (f"Expected {stage} duration as a proportion of SPT"
                if variable == "HD_DENS" else f"Expected {stage} duration from posterior mass (minutes)")
    if variable in {"HD_DENS", "HD_MINS"} and "MIX" in factors:
        state_set = factors["MIX"].replace("_", "/")
        return (f"Expected time jointly ambiguous across {state_set}, as a proportion of SPT"
                if variable == "HD_DENS" else f"Expected minutes jointly ambiguous across {state_set}")
    base = METRIC_BRIEFS.get(variable, variable.replace("_", " ").title())
    context = stratum_brief(stratum)
    return base if not context else f"{base} ({context})"


def feature_label(variable: str, stratum: str) -> str:
    """Return a compact display label, independent of the panel using it."""
    base = METRIC_LABELS.get(variable, METRIC_BRIEFS.get(variable, variable.replace("_", " ").title()))
    if not stratum or stratum == ".":
        return base
    contexts: list[str] = []
    variation = ""
    window = ""
    transition = ""
    mixture = ""
    stage = ""
    ultradian_quartile = ""
    ultradian_metric = ""
    for field in stratum.split(";"):
        factor, level = field.split("/", 1)
        if factor == "SS":
            stage = level
            contexts.append(level)
        elif factor == "MIX":
            mixture = level.replace("_", "/")
        elif factor in {"TRANS", "BOUT_TRANS"}:
            transition = level.replace("to", "->")
            contexts.append(transition)
        elif factor == "STABLE":
            contexts.append("stable" if level == "1" else "unstable")
        elif factor == "LRLAG":
            contexts.append(duration_label(level))
        elif factor == "LRWIN":
            window = duration_label(level)
            contexts.append(f"window={window}")
        elif factor == "VAR":
            variation = level.replace("_", "/")
        elif factor == "UQ":
            ultradian_quartile = level
        elif factor == "UMETRIC":
            ultradian_metric = level
        elif factor == "UVAR":
            variation = level.replace("_", "/")
        # NSS is implicit in the associated staging/metric convention and is
        # deliberately omitted from compact viewer labels.
    lower_after_prefix = variable in {
        "SLOPE", "U_MOD", "U_MOD_COS", "U_MOD_SIN", "U_PARTIAL_R2", "U_R2", "U_SLOPE",
    }
    if ultradian_metric:
        return f"{ultradian_metric} {base}" + (f" ({ultradian_quartile})" if ultradian_quartile else "")
    if variation:
        if variation == "H":
            variation = "Entropy"
        text = base[:1].lower() + base[1:] if lower_after_prefix else base
        if variable in {"BETWEEN_SD", "WITHIN_SD", "SLOPE"} and variation in {"W", "N1", "N2", "N3", "R", "NR"}:
            variation = f"P({variation})"
        return f"{variation} {text}" + (f" ({window})" if window else "")
    if mixture:
        return f"{mixture} {base}"
    if variable in {"HD_DENS", "HD_MINS", "HD_PCT"} and stage:
        if variable == "HD_DENS":
            return f"{stage} density"
        if variable == "HD_MINS":
            return f"{stage} minutes"
        return f"{stage} %"
    if transition:
        return f"{transition} {base[:1].lower()}{base[1:]}"
    return base if not contexts else f"{base}({', '.join(contexts)})"


def label_qualifier(stratum: str) -> str:
    """Return the sole allowed stage qualifier for a duplicate label."""
    for field in stratum.split(";"):
        factor, level = field.split("/", 1)
        if factor == "NSS" and level == "3":
            return "3-stage"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path, nargs="?",
                        help="generated C++ header (omit with --json-output)")
    parser.add_argument("--json-output", type=Path,
                        help="also write the equivalent self-contained V1 JSON layout")
    args = parser.parse_args()
    if not args.output and not args.json_output:
        parser.error("provide an output header and/or --json-output")

    features: list[tuple[str, str, str]] = []
    components: list[tuple[str, str, int, int, str]] = []
    seen_features: set[tuple[str, str, str]] = set()
    seen_components: set[str] = set()

    for line_no, raw in enumerate(args.source.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#") or line == "DOMAIN\tSTRATUM\tVARIABLE":
            continue
        fields = raw.split("\t")
        if fields[0] == "COMPONENT":
            if len(fields) != 6:
                raise SystemExit(f"{args.source}:{line_no}: COMPONENT needs 6 columns")
            ident, panel, pc, direction, label = fields[1:]
            try:
                pc_number, sign = int(pc), int(direction)
            except ValueError as error:
                raise SystemExit(f"{args.source}:{line_no}: invalid component number/sign") from error
            if not ident or not panel or pc_number < 1 or sign not in {-1, 1} or not label:
                raise SystemExit(f"{args.source}:{line_no}: invalid COMPONENT")
            if ident in seen_components:
                raise SystemExit(f"{args.source}:{line_no}: duplicate COMPONENT {ident}")
            seen_components.add(ident)
            components.append((ident, panel, pc_number, sign, label))
        else:
            if len(fields) != 3:
                raise SystemExit(f"{args.source}:{line_no}: feature needs 3 columns")
            feature = tuple(fields)
            if feature in seen_features:
                raise SystemExit(f"{args.source}:{line_no}: duplicate feature {feature}")
            seen_features.add(feature)
            features.append(feature)

    if not features or not components:
        raise SystemExit(f"{args.source}: features and COMPONENTs are required")

    lines = [
        "// Generated by scripts/generate_hdstats_pc_schema.py; do not edit.",
        "#ifndef __HDSTATS_PC_SCHEMA_H__",
        "#define __HDSTATS_PC_SCHEMA_H__",
        "namespace hdstats_pc_generated {",
        "static const generated_feature_t features[] = {",
    ]
    lines += [f"  {{ {cpp(panel)}, {cpp(stratum)}, {cpp(variable)} }}," for panel, stratum, variable in features]
    lines += [
        "};",
        "static const int n_features = sizeof(features) / sizeof(features[0]);",
        "static const generated_component_t components[] = {",
    ]
    lines += [
        f"  {{ {cpp(ident)}, {cpp(panel)}, {pc}, {sign}, {cpp(label)} }},"
        for ident, panel, pc, sign, label in components
    ]
    lines += [
        "};",
        "static const int n_components = sizeof(components) / sizeof(components[0]);",
    ]
    lines += [
        "}",
        "#endif",
        "",
    ]
    if args.output:
        args.output.write_text("\n".join(lines))

    if args.json_output:
        panels: dict[str, dict] = {}
        panel_codes: dict[str, str] = {}
        code_panels: dict[str, str] = {}
        for ident, panel, _pc, _sign, _label in components:
            match = re.fullmatch(r"PC_([A-Z][A-Z0-9]*)_[1-9][0-9]*", ident)
            if not match:
                raise SystemExit(f"component ID must encode a panel code (PC_CODE_N): {ident}")
            code = match.group(1)
            if panel in panel_codes and panel_codes[panel] != code:
                raise SystemExit(f"panel {panel} has conflicting component codes")
            if code in code_panels and code_panels[code] != panel:
                raise SystemExit(f"component code {code} is shared by multiple panels")
            panel_codes[panel] = code
            code_panels[code] = panel
        descriptions: dict[str, dict[str, str]] = {}
        description_strata: dict[str, str] = {}
        for panel, stratum, variable in features:
            if panel not in panel_codes:
                raise SystemExit(f"panel {panel} is missing a coded COMPONENT entry")
            code = panel_codes[panel]
            panels.setdefault(code, {"id": code, "features": [], "n_components": 0})
            identifier = feature_key(variable, stratum)
            panels[code]["features"].append(identifier)
            descriptions.setdefault(identifier, {
                "label": feature_label(variable, stratum),
                "description": feature_brief(variable, stratum),
            })
            description_strata.setdefault(identifier, stratum)
        by_label: dict[str, list[str]] = {}
        for identifier, metadata in descriptions.items():
            by_label.setdefault(metadata["label"], []).append(identifier)
        for identifiers in by_label.values():
            if len(identifiers) > 1:
                for identifier in identifiers:
                    qualifier = label_qualifier(description_strata[identifier])
                    if qualifier:
                        descriptions[identifier]["label"] += f" [{qualifier}]"
        for metadata in descriptions.values():
            label = metadata["label"]
            metadata["label"] = label[:1].upper() + label[1:]
            if not metadata["label"][0].isalpha() or not metadata["label"][0].isupper():
                raise SystemExit(f"feature label must start with a capital letter: {metadata['label']}")
        labels = [metadata["label"] for metadata in descriptions.values()]
        if len(labels) != len(set(labels)):
            raise SystemExit("generated HDSTATS feature labels are not unique")
        for ident, panel, pc, sign, label in components:
            code = panel_codes[panel]
            panels[code]["n_components"] += 1
        document = {
            "format": "luna-hdstats-pc-spec",
            "version": 1,
            "feature_universe": "HDSTATS-PC-FEATURES-V1",
            "component_id_style": "r_safe",
            "feature_descriptions": [
                {"id": identifier, **metadata}
                for identifier, metadata in sorted(descriptions.items())
            ],
            "panels": list(panels.values()),
        }
        args.json_output.write_text(json.dumps(document, indent=2) + "\n")


if __name__ == "__main__":
    main()
