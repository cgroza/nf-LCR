#!/usr/bin/env bash
set -euo pipefail

image_tag="${1:-cgroza/nf-chr22q:wfmash-0.24.2}"
ca_bundle="$(mktemp -t nf-chr22q-ca)"
trap 'rm -f "$ca_bundle"' EXIT

command -v docker >/dev/null || {
  printf 'Docker is required.\n' >&2
  exit 1
}
command -v security >/dev/null || {
  printf 'This build script requires the macOS security command.\n' >&2
  exit 1
}

security find-certificate -a -p \
  /System/Library/Keychains/SystemRootCertificates.keychain > "$ca_bundle"
security find-certificate -a -p \
  /Library/Keychains/System.keychain >> "$ca_bundle"

docker buildx build \
  --platform linux/amd64 \
  --secret "id=mac_ca,src=$ca_bundle" \
  --load \
  --tag "$image_tag" \
  .

docker image inspect "$image_tag" --format '{{.Os}}/{{.Architecture}} {{.Id}}'
