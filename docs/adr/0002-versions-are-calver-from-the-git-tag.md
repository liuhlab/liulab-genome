# 2. Versions are CalVer `YYYY.MM.MICRO`, derived from the git tag

hatch-vcs derives the version from the tag, so a tag *is* a release and a version is never
hand-edited — there is no second copy of the number to drift out of agreement with the first. What it
costs: no semantic signal for breakage, so a consumer cannot read compatibility off the number and
has to read `CHANGELOG.md`.
