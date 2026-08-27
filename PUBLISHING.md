# Private/public versioning

## Decision

Use two repositories in two directories, each with its own `origin` and an
independent Git history.  Do not use private and public branches in one
repository.

A Git branch is not a confidentiality boundary.  Objects reachable only
from a private branch can still be disclosed by pushing the wrong branch or
tag, and rewriting the public branch later does not reliably retract them.
The private and public repositories also have different provenance,
licensing, and disclosure requirements, so independent histories are the
cleaner model.

The private repository is canonical.  The public repository is a deliberate
projection defined by `PUBLIC_FILES.txt`; it is not a mirror of the complete
private tree.  The export command copies only allowlisted paths, refuses
secret-like or Git-internal paths, records a managed-file manifest, and
writes `CHECKSUMS.sha256` for the exported payload.

## One-time setup

1. Keep the development working tree as the private repository and give it
   a private `origin`.
2. Create a separate sibling directory, initialize a new Git repository, and
   give it the public `origin`.
3. Do not copy `.git`, tags, bundles, patches, reflogs, or private commit IDs
   into the public repository.

## Release checklist

1. Commit and verify the private state.
2. Run `make verify-all` and `make paper`.
3. Export into the independent public working tree:

   ```text
   make export-public EXPORT_DIR=/absolute/path/to/public-repository
   ```

4. In the public repository, inspect `git status`, `git diff`,
   `.public-export-manifest.json`, and `CHECKSUMS.sha256`.
5. Re-run `make verify-all` and `make paper` in the public repository.
6. Commit and push from the public repository only after that review.

The exporter treats existing managed public files as generated outputs.  It
refuses to overwrite a managed file that has been edited since the previous
export.  Public-only, unmanaged files are preserved.  A file removed from
the allowlist is removed on the next export only if it still matches the
previously exported checksum.

## Remotes

Give each repository only the `origin` it is meant to push to.  If a
read-only public remote is useful in the private repository, configure it
with an intentionally invalid push URL.  This is optional; the normal
release path should still be the explicit export into the public directory.
