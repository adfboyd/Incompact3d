#!/bin/bash
# Recursively touch all files in a directory tree to reset their access/mod times.
# Useful on NVMe scratch filesystems that auto-delete old files.
#
# Usage: ./touch_all.sh [directory]
# Default: current directory

TARGET="${1:-.}"

if [[ ! -d "$TARGET" ]]; then
    echo "Error: '$TARGET' is not a directory" >&2
    exit 1
fi

total=$(find "$TARGET" -type f | wc -l)

if [[ "$total" -eq 0 ]]; then
    echo "No files found in '$TARGET'" >&2
    exit 0
fi

bar_width=40
done=0
batch=200
start_time=$SECONDS

fmt_time() {
    local s=$1
    printf '%d:%02d:%02d' $(( s/3600 )) $(( (s%3600)/60 )) $(( s%60 ))
}

draw_bar() {
    local pct=$(( done * 100 / total ))
    local filled=$(( done * bar_width / total ))
    local bar=$(printf '%*s' "$filled" '' | tr ' ' '#')
    local empty=$(printf '%*s' $(( bar_width - filled )) '')
    local elapsed=$(( SECONDS - start_time ))
    local eta_str
    if [[ $done -gt 0 && $done -lt $total ]]; then
        local eta=$(( (total - done) * elapsed / done ))
        eta_str="ETA $(fmt_time $eta)"
    elif [[ $done -ge $total ]]; then
        eta_str="done in $(fmt_time $elapsed)"
    else
        eta_str="ETA --:--:--"
    fi
    printf "\r[%s%s] %3d%% (%d/%d)  elapsed %s  %s" \
        "$bar" "$empty" "$pct" "$done" "$total" \
        "$(fmt_time $elapsed)" "$eta_str"
}

find "$TARGET" -type f -print0 \
  | xargs -0 -n "$batch" sh -c '
        touch "$@"
        echo ${#@}
    ' _ \
  | while IFS= read -r n; do
        (( done += n ))
        draw_bar
    done

done=$total
draw_bar
echo
echo "Done — touched $total files in '$TARGET'" >&2
