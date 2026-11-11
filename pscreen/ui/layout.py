"""Viewport-aware sizing: a zero-height component reports the browser window height so pages can size panels to fit without page scrolling."""
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components

_vh = components.declare_component("pscreen_vh", path=str(Path(__file__).parent / "vh_component"))


def report_viewport():
    """Call once per run, near the top of the app."""
    with st.sidebar:                                  # the sidebar is hidden, so this adds no visible gap
        v = _vh(key="vh", default=None)
    if v: st.session_state["_vh"] = int(v)


def viewport_h(default: int = 860) -> int:
    return int(st.session_state.get("_vh", default))
