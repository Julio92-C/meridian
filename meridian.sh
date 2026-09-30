#!/usr/bin/env bash
#
# meridian.sh — run the MERIDIAN container like the native pipeline:
#
#     ./meridian.sh <config.yaml>
#
# Wraps the container so you only pass a config, exactly as with the native
#
#     Rscript run_pipeline.R <config.yaml>
#
# The config's `project_root:` — an ABSOLUTE host path holding Metadata/ and
# Reports/ — is bound into the container at the same path, so the config needs
# no edits and outputs land under <project_root>/test_run/.
#
# Runtime is auto-detected: Apptainer/Singularity if present (HPC), else Docker
# (laptop). Overrides via environment:
#   MERIDIAN_SIF     path to the .sif        (default: <script dir>/meridian_v1.0.0.sif)
#   MERIDIAN_IMAGE   docker/oras image ref   (default: ghcr.io/julio92-c/meridian:v1.0.0)
#
set -euo pipefail

IMAGE_DEFAULT="ghcr.io/julio92-c/meridian:v1.0.0"

# --- args ---------------------------------------------------------------------
CONFIG="${1:-}"
[ -n "$CONFIG" ]  || { echo "usage: meridian.sh <config.yaml>" >&2; exit 2; }
[ -f "$CONFIG" ]  || { echo "meridian.sh: config not found: $CONFIG" >&2; exit 2; }
CONFIG="$(readlink -f "$CONFIG")"
CONFDIR="$(dirname "$CONFIG")"
HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

# --- resolve project_root from the config (must be an absolute host path) ------
PROOT="$(grep -E '^[[:space:]]*project_root:' "$CONFIG" | head -1 \
         | sed -E 's/[^:]*:[[:space:]]*"?([^"#]+)"?.*/\1/' | sed -E 's/[[:space:]]+$//')"
[ -n "$PROOT" ] || { echo "meridian.sh: no project_root found in $CONFIG" >&2; exit 1; }
case "$PROOT" in
  /*) : ;;
  *)  echo "meridian.sh: project_root must be an ABSOLUTE host path for the container (got: '$PROOT')." >&2
      echo "  Edit $CONFIG so project_root points at the directory holding Metadata/ and Reports/." >&2
      exit 1 ;;
esac
[ -d "$PROOT" ] || { echo "meridian.sh: project_root does not exist: $PROOT" >&2; exit 1; }

# Bind project_root; add the config dir too, unless it already lives under it.
BINDS=("$PROOT:$PROOT")
case "$CONFDIR/" in
  "$PROOT"/*) : ;;
  *) BINDS+=("$CONFDIR:$CONFDIR") ;;
esac

# --- pick a runtime and exec --------------------------------------------------
run_apptainer() {   # $1 = apptainer|singularity
  local rt="$1"
  local sif="${MERIDIAN_SIF:-$HERE/meridian_v1.0.0.sif}"
  if [ ! -f "$sif" ]; then
    echo "meridian.sh: image not found: $sif" >&2
    echo "  build it once with:  $rt build \"$sif\" docker://${MERIDIAN_IMAGE:-$IMAGE_DEFAULT}" >&2
    echo "  or point MERIDIAN_SIF at an existing .sif." >&2
    exit 1
  fi
  # The report stage renders templates/report.qmd IN PLACE (Quarto writes its
  # intermediate .md / _files next to the .qmd), and the qmd locates the repo
  # via rprojroot -> it must live inside /opt/meridian. But the image rootfs is
  # read-only. So stage a writable copy of templates/ and bind it over the
  # image's templates dir; the rest of the repo stays read-only and readable.
  local tmpl; tmpl="$(mktemp -d "${TMPDIR:-/tmp}/meridian_tmpl.XXXXXX")"
  trap 'rm -rf "$tmpl"' EXIT
  if ! "$rt" exec "$sif" cp -a /opt/meridian/templates/. "$tmpl/" 2>/dev/null; then
    echo "meridian.sh: warning: could not stage a writable templates/ dir; the HTML report may fail on a read-only image." >&2
  fi
  local args=(run --pwd /opt/meridian -B "$tmpl:/opt/meridian/templates")
  for b in "${BINDS[@]}"; do args+=(-B "$b"); done
  echo ">> $rt: $sif" >&2
  "$rt" "${args[@]}" "$sif" "$CONFIG"   # not exec, so the EXIT trap cleans up
  exit $?
}

run_docker() {
  local img="${MERIDIAN_IMAGE:-$IMAGE_DEFAULT}"
  local args=(--rm -w /opt/meridian)
  for b in "${BINDS[@]}"; do args+=(-v "$b"); done
  echo ">> docker: $img" >&2
  exec docker run "${args[@]}" "$img" "$CONFIG"
}

if   command -v apptainer   >/dev/null 2>&1; then run_apptainer apptainer
elif command -v singularity >/dev/null 2>&1; then run_apptainer singularity
elif command -v docker      >/dev/null 2>&1; then run_docker
else
  echo "meridian.sh: no container runtime found (need apptainer, singularity, or docker)." >&2
  exit 1
fi
