"""Exercise all tool calls over a real MCP stdio connection in a fresh process."""

from __future__ import annotations

import asyncio
import json
import sys

import numpy as np
import pytest

from spatial_audit.schema import write_json


def test_stdio_all_tools_and_error_recovery(tmp_path):
    mcp = pytest.importorskip("mcp")
    from mcp.client.stdio import stdio_client

    xy = np.indices((4, 4)).reshape(2, -1).T
    files = []
    for i in range(4):
        name = f"p{i}.npz"
        np.savez(
            tmp_path / name,
            X=np.arange(16.0)[:, None] + i,
            coords=xy * 100.0,
            features=["intensity"],
            C=xy[:, [0]],
        )
        files.append(
            {
                "path": name,
                "patient_id": f"person-{i}",
                "group": str(i // 2),
                "slide_id": str(i % 2),
                "batch": "shared",
                "stage": [i // 2, 0.0],
            }
        )
    manifest = tmp_path / "manifest.json"
    write_json(
        manifest,
        {
            "cohorts": [
                {"id": "synthetic-protocol", "kind": "synthetic", "units": "um", "files": files}
            ]
        },
    )

    async def check():
        params = mcp.StdioServerParameters(
            command=sys.executable, args=["-m", "spatial_audit.mcp_server"], env={"PYTHONUTF8": "1"}
        )
        async with stdio_client(params) as streams, mcp.ClientSession(*streams) as sess:
            await sess.initialize()
            listed = await sess.list_tools()
            assert {x.name for x in listed.tools} == {
                "compute_exact_p_floor",
                "simulate_trap_cohort",
                "audit_cohort",
                "check_design_collinearity",
            }

            async def run(name, args):
                response = await sess.call_tool(name, args)
                assert not response.isError
                contents = response.structuredContent
                assert isinstance(contents, dict)
                json.dumps(contents, allow_nan=False)
                return contents

            args = {"manifest": str(manifest), "cohort_id": "synthetic-protocol"}
            audit = await run("audit_cohort", args)
            assert audit["n_patients"] == 4 and len(audit["dimensions"]) == 5
            result = await run("check_design_collinearity", args)
            assert result["code"] == "UNIDENTIFIABLE_DESIGN"
            failed = await sess.call_tool("compute_exact_p_floor", {"n0": 0, "n1": 3})
            assert failed.isError
            floor = await run("compute_exact_p_floor", {"n0": 3, "n1": 3})
            assert floor["p_floor"] == 0.1
            sim = await run("simulate_trap_cohort", {"trap": "stage"})
            assert sim["kind"] == "synthetic" and sim["dimensions"]["design"]["status"] == "red"

    asyncio.run(asyncio.wait_for(check(), timeout=60))
