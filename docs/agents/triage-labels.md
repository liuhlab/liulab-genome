# Triage labels

The skills speak of five triage roles, and this tracker carries exactly those five. Each
label string equals its role name, so the mapping is the identity — there is no
translation step to get wrong, and nothing here to keep in sync with a second column.

| Label | Meaning |
| --- | --- |
| `needs-triage` | Maintainer needs to evaluate this issue |
| `needs-info` | Waiting on reporter for more information |
| `ready-for-agent` | Spec or ticket is agent-ready |
| `ready-for-human` | Needs a human: judgment call, design decision, or external access |
| `wontfix` | This will not be worked on |

Apply one with `gh issue edit <number> --add-label "ready-for-agent"`.

Triage is for issues someone else wrote. Spec output — the tickets `/to-spec` emits — is
agent-ready as written and skips triage.
