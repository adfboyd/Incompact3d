#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE' >&2
Usage: new_run.sh RUN_NAME INPUT

Create ellipsoid_dev/runs/RUN_NAME, copy INPUT to input.i3d there, and write a
local MPI run script.

INPUT may be either:
  - a filename from ellipsoid_dev/inputs, e.g. input_sphere_inviscid_rotating_test.i3d
  - an explicit path to an input file

Use NP to override the local MPI rank count:
  NP=4 ./run_local.sh
USAGE
}

if [[ $# -ne 2 ]]; then
    usage
    exit 2
fi

run_name=$1
input_arg=$2

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
dev_dir=$(cd "$script_dir/.." && pwd)
repo_root=$(cd "$dev_dir/.." && pwd)
run_dir="$dev_dir/runs/$run_name"

if [[ "$run_name" == *"/"* || "$run_name" == "." || "$run_name" == ".." ]]; then
    echo "Error: RUN_NAME must be a simple directory name" >&2
    exit 2
fi

if [[ -f "$input_arg" ]]; then
    input_path=$(cd "$(dirname "$input_arg")" && pwd)/$(basename "$input_arg")
elif [[ -f "$dev_dir/inputs/$input_arg" ]]; then
    input_path="$dev_dir/inputs/$input_arg"
else
    echo "Error: input not found: $input_arg" >&2
    echo "Looked in current path and $dev_dir/inputs" >&2
    exit 1
fi

mkdir -p "$run_dir"
cp "$input_path" "$run_dir/input.i3d"

cat > "$run_dir/run_local.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail

NP=\${NP:-8}
mpirun --oversubscribe -np "\$NP" "$repo_root/build/bin/xcompact3d" input.i3d
EOF
chmod +x "$run_dir/run_local.sh"

cat > "$run_dir/README.md" <<EOF
# $run_name

Input copied from:

\`\`\`
$input_path
\`\`\`

Run locally from this directory:

\`\`\`sh
./run_local.sh
\`\`\`

Override MPI ranks with:

\`\`\`sh
NP=4 ./run_local.sh
\`\`\`
EOF

echo "Created $run_dir"
