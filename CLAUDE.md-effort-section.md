# Effort & usage management (Pro plan)

I'm on a Claude Pro plan, which has the tightest rolling usage limits of any
paid tier. Manage your own effort accordingly — the goal is to get real work
done without burning through my 5-hour/weekly budget on wasted tool calls,
redundant context, or oversized reasoning for small tasks.

## Subagents (Task tool)

Delegate to a subagent when a piece of work can be done in isolation and only
the *result* matters to the main thread — codebase exploration, searching for
a symbol/pattern across many files, running and summarizing a test suite,
independent verification/review of a change, or research that would otherwise
dump a lot of raw text into my main context.

Keep work inline (don't spawn a subagent) when steps are tightly coupled,
require back-and-forth judgment calls, or the task is already small enough
that spawning a subagent (fresh context load: this file, git status, skills)
would cost more than it saves.

Default subagent model: **Sonnet**, not Haiku, at **medium effort**, for
everything except the smallest, most mechanical lookups (e.g. "does this repo
use pnpm or npm" — trivial enough for low effort). Reserve high/max effort for
subagents doing genuinely hard reasoning (tricky debugging, architectural
tradeoffs), and only when the task clearly warrants it.

Don't fan out many subagents in parallel by default — each one is a full
fresh context load and a real chunk of budget. Parallelize only when the
subtasks are genuinely independent and the task is large enough that serial
execution would clearly take too long; otherwise prefer one or two targeted
subagents over a swarm. If you think a task would benefit from heavier
parallel/swarm-style delegation than this default, say so and ask before
doing it.

## Context and effort hygiene

- Prefer `/compact` proactively once a session's context is getting large,
  rather than letting it balloon — pass instructions on what to keep (e.g.
  "keep the diff and open questions, drop tool output").
- Keep this file lean. If you find yourself wanting to add a lot of detailed,
  situational instructions here, propose a skill instead — skills load on
  demand, this file loads on every session start.
- When a tool would return a lot of verbose, low-signal output (full test
  logs, large file dumps), summarize or filter it down rather than echoing it
  all into context, and prefer a hook for this where one exists.

## Usage limits — what you can and can't check yourself

You (Claude Code) cannot query my remaining Pro plan quota programmatically —
there's no API or file you can read for this. The only view into it is the
`/usage` slash command, which is interactive and shows the current session's
token/cost breakdown plus my plan's usage bars, but only when I run it or
paste you the output.

So:
- If a session has been running long, or you're about to kick off something
  expensive (a big parallel delegation, a large refactor, heavy extended
  thinking), say so and suggest I run `/usage` to check headroom before
  proceeding — don't assume you know where I stand.
- If I paste you `/usage` output, use it to calibrate: back off toward
  low-effort/serial/inline work if I'm close to a limit, and don't ask me to
  check again immediately after.
- Never claim to know my remaining usage unless I've just given it to you.
