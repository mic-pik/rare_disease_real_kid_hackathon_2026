"""Parse local model tool calls and user-channel answers."""

import re

CALL_RE = re.compile(
    r'<atem:invoke\s+name="(?P<name>[^"]+)">(?P<body>.*?)</atem:invoke>', re.S)

PARAM_RE = re.compile(
    r'<atem:parameter\s+name="(?P<key>[^"]+)">(?P<val>.*?)</atem:parameter>', re.S)

def parse_calls(text: str) -> list[dict]:
    out = []
    for m in CALL_RE.finditer(text):
        args = {}
        for p in PARAM_RE.finditer(m.group("body")):
            v = p.group("val").strip()
            if v.lower() in ("true", "false"):
                v = v.lower() == "true"
            elif re.fullmatch(r"-?\d+", v):
                v = int(v)
            args[p.group("key")] = v
        out.append({"name": m.group("name"), "arguments": args})
    return out

ANY_CHANNEL = re.compile(r"(?:assistant\s*)?to=(?P<chan>user|self)\s*", re.I)

def user_answer(text: str) -> str | None:
    """The user-facing part of a generation, or None if it never got there.

    Returns every `to=user` segment joined in order, with `to=self` reasoning
    dropped. None means the generation never reached the user channel at all,
    which the caller treats as "not finished yet" rather than as an answer.
    """
    text = text or ""
    marks = list(ANY_CHANNEL.finditer(text))
    if not marks:
        return None
    parts = []
    for i, m in enumerate(marks):
        if m.group("chan").lower() != "user":
            continue
        end = marks[i + 1].start() if i + 1 < len(marks) else len(text)
        seg = text[m.end():end].strip()
        if seg:
            parts.append(seg)
    if not parts:
        # Reached the user channel but said nothing there yet.
        return None
    return "\n\n".join(parts)
